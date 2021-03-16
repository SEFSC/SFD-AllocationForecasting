remotes::install_github(repo="r4ss/r4ss",ref="development")

library(r4ss)

#This is currently set up for yellowtail snapper. Will need work to generalize to all assessment formats.

run.projections<-function(assessment_dir, #Here you set the location of a previously fit SS stock assessment
                          Target.Depletion=0.3, #Set your target OFL depletion proxy 
                          PstarFraction=0.75, #Set the P* target as a fraction of the OFL target
                          Recent_Catch=matrix(data=c(898.0704387,995.201,113.28,120.76,1696.550856,805.737),nrow=2,ncol=3,byrow=FALSE), #Recent Catches for 2018 and 2019
                          Allocate_Catch=FALSE, Comm.Frac=NULL, Rec.Frac=NULL, #Should catch have a fixed allocation between rec and commercial and if so what are those ratios.
                          Catch.Treshold=0.00001, Depletion.Threshold=0.0001, Annual.F.Threshold=0.0001, #These are all just thresholds for when to stop adjusting catches
                          Allocation.Threshold=0.0001, Relative.F.Threshold=0.0001) #increase them if too slow to converge.
{
  #Set large max print to avoid issues with writing out a large forecast file.
  options(max.print=1000000)
  #First set up a working director for running projections in (to avoid overwriting the base files with a failed model run)
  #then copy all of the assessment files to this working folder (ignore any output directories that have been previously created)
  if(dir.exists(paste0(assessment_dir,"/Working_dir"))){
    unlink(paste0(assessment_dir,"/Working_dir"),recursive = TRUE)
  }
  dir.create(paste0(assessment_dir,"/Working_dir"))
  temp.files<-list.files(path=assessment_dir)
  folders<-c(which(temp.files=="Depletion_target"),which(temp.files=="Pstar_target"),which(temp.files=="Working_dir"))
  if(length(folders)>0){
    temp.files<-temp.files[-folders]
  }
  file.copy(from=paste0(assessment_dir,'/',temp.files),to=paste0(assessment_dir,"/Working_dir/",temp.files))
  
  #Set working directory and read in the assessment files
  setwd(paste0(assessment_dir,"/Working_dir"))
  start<-SS_readstarter()
  dat<-SS_readdat(file=start$datfile,version=3.3)
  ctl<-SS_readctl(file=start$ctlfile,version=3.3,use_datlist = TRUE, datlist = dat)
  results<-SS_output(dir=getwd(),covar=FALSE,checkcor = FALSE)
  forecast<-SS_readforecast() 
  parlist<-SS_readpar_3.30(parfile="ss.par",datsource = dat, ctlsource = ctl)
  
  #Modify assessment files to produce expected timeseries length and outputs
  #100 year projections allow equilibrium to be achieved
  forecast$Nforecastyrs<-100
  parlist$recdev_forecast<-matrix(NA,nrow=100,ncol=2)
  parlist$Fcast_impl_error<-matrix(NA,nrow=100,ncol=2)
  #Need to set recdevs and implementation error to 0 for all projection years 
  #so that reading from par file is possible
  parlist$recdev_forecast[,1]<-2018:2117
  parlist$Fcast_impl_error[,1]<-2018:2117
  parlist$recdev_forecast[,2]<-rep(0,100)
  parlist$Fcast_impl_error[,2]<-rep(0,100)
  colnames(parlist$recdev_forecast)<-c("year","recdev")
  colnames(parlist$Fcast_impl_error)<-c("year","impl_error")
  
  #Adjust the starter file to read from par file, perform no fitting (This should already have been done),
  #and set the depletion value to be relative to unexploited biomass and have no scaling 
  #(so that correct depletion target can be found).
  start$init_values_src<-1
  start$last_estimation_phase<-0
  start$depl_basis<-1
  start$depl_denom_frac<-1
  
  
  #Get the timeseries of historic/projected catches and F
  TimeSeries<-results$timeseries
  #Set all future projections to fish at constant apical F that matches recent years
  #get the years of timeseries F's based no the forecast year range
  TargetYears<-TimeSeries[TimeSeries$Yr>=forecast$Fcast_years[3] & TimeSeries$Yr<=forecast$Fcast_years[4],]
  RecentYears<-TimeSeries[TimeSeries$Yr>=2018 & TimeSeries$Yr<=2019,]
  Recent_Catch_Achieved<-RecentYears[,c(14,25,33)]
  Starting_Adjust<-Recent_Catch/Recent_Catch_Achieved
  F_by_Fleet<-apply(TargetYears[,c(2,19+0:2*8)],2,mean)
  
  forecast$InputBasis<-99
  forecast_F<-matrix(1,nrow=300,ncol=4)
  forecast_F[,1]<-sort(rep(2018:2117,3))
  forecast_F[,3]<-rep(1:3,100)
  forecast_F[1:6,4]<-c(unlist(RecentYears[1,19+(0:2)*8]),unlist(RecentYears[2,19+(0:2)*8]))*c(unlist(Starting_Adjust[1,]),unlist(Starting_Adjust[2,]))
  forecast_F[-c(1:6),4]<-rep(F_by_Fleet[2:4],98)
  forecast_F<-as.data.frame(forecast_F)
  colnames(forecast_F)<-c("Year","Seas","Fleet","Catch or F")
  forecast$ForeCatch<-forecast_F
  #Save all the modified files and then perform a base run of SS so that output is dimensioned correctly with 
  #a 100 year projection series
  SS_writepar_3.30(parlist = parlist,outfile="ss.par",overwrite = TRUE)
  SS_writeforecast(mylist=forecast,overwrite = TRUE)
  SS_writestarter(mylist=start,overwrite = TRUE)
  
  shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
  
  #Now start a loop of projecting and modifying fixed F's until the desired 
  #catch projections are achieved
  keepFitting<-TRUE
  fitting_depletion<-TRUE
  while(keepFitting){
    #Read in the SS results for landings and stock status to determine if desired
    #targets have been achieved
    resultsFit<-SS_output(dir=getwd(),covar=FALSE,checkcor = FALSE)
    TimeFit<-resultsFit$timeseries
    SPRfit<-resultsFit$sprseries
    plot(SPRfit$F_report[SPRfit$Yr>=2020])
    Recent_Catch_Achieved<-TimeFit[TimeFit$Yr>=2018 & TimeFit$Yr<=2019,c(14,25,33)]
    
    if(fitting_depletion==TRUE){
      #Check the equilibrium depletion relative to the desired target and develop
      #an adjustment factor based on this to scale all F's in next projection loop.
      #Add randomization to the value to stop the model from getting caught in local
      #minimums or oscillating to speed convergence. 
      Depletion<-SPRfit$Deplete
      DepletionScale<-runif(1,0,0.25)*((mean(Depletion[(length(Depletion)-10):length(Depletion)])/Target.Depletion)-1)+1
      if(abs(1-DepletionScale)>=Depletion.Threshold){keepFitting<-TRUE}else{keepFitting<-FALSE}
      #Calculate the average F at equilibrium that all F's will be scaled to in order
      #to achieve equal F in every year. As depletion approaches the target value this 
      #F will approach F(OFL).
      F_report<-SPRfit$F_report
      FScale<-mean(F_report[(length(F_report)-10):length(F_report)])
      Fpstar<-FScale
    }else{
      keepFitting<-FALSE
      DepletionScale<-1
      FScale<-PstarFraction*Fpstar
    }
    
    for(i in 1:2){
      for(j in 1:3){
        Catch_scale<-runif(1,0,0.5)*(Recent_Catch[i,j]/Recent_Catch_Achieved[i,j]-1)+1
        forecast$ForeCatch[((i-1)*3+j),4]<-forecast$ForeCatch[((i-1)*3+j),4]*Catch_scale
        if(abs(1-Catch_scale)<=Catch.Treshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
      }
    }
    #Now loop over each year of the projection and calculate adjustments based on allocation fraction and annual F.
    for(i in 2020:2117){
      #Select the forecast catch row by excluding the first 2 years of recent catch 
      row<-6+3*(i-2020)
      
      #Calculate the annual F multiplier by comparing this years value to the equilibrium value and adjust to match.
      Fmult1<-runif(1,0,0.5)*(FScale/SPRfit$F_report[SPRfit$Yr==i]-1)+1
      if(abs(1-Fmult1)<=Annual.F.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
      
      #If allocating catch calculate the scaling needed to bring allocations to expected value
      if(Allocate_Catch){
        group1Catch<-TimeFit$`retain(B):_1`[TimeFit$Yr==i]
        group2Catch<-sum(TimeFit$`retain(B):_2`[TimeFit$Yr==i],TimeFit$`retain(B):_3`[TimeFit$Yr==i])
        rand1<-runif(1,0,0.5)
        Fmult2G1<-rand1*(Comm.Frac/(group1Catch/(group1Catch+group2Catch))-1)+1
        if(abs(1-Fmult2G1)<=Allocation.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        Fmult2G2<-rand1*(Rec.Frac/(group2Catch/(group1Catch+group2Catch))-1)+1
        if(abs(1-Fmult2G2)<=Allocation.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        
        #Scale the relative F's of fleets if they are not under a fixed allocation scheme. 
        rand2<-runif(1,0,0.5)
        
        Fmult3F1<-1
        Fmult3F2<-rand2*((F_by_Fleet[3]/(F_by_Fleet[3]+F_by_Fleet[4]))/(TimeFit$`F:_2`[TimeFit$Yr==i]/(TimeFit$`F:_2`[TimeFit$Yr==i]+TimeFit$`F:_3`[TimeFit$Yr==i]))-1)+1
        if(abs(1-Fmult3F2)<=Relative.F.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        Fmult3F3<-rand2*((F_by_Fleet[4]/(F_by_Fleet[3]+F_by_Fleet[4]))/(TimeFit$`F:_3`[TimeFit$Yr==i]/(TimeFit$`F:_2`[TimeFit$Yr==i]+TimeFit$`F:_3`[TimeFit$Yr==i]))-1)+1
        if(abs(1-Fmult3F3)<=Relative.F.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        
      }else{
        Fmult2G1<-1
        Fmult2G2<-1
        #Scale the relative F's of all fleets.
        rand2<-runif(1,0,0.5)
        Fmult3F1<-rand2*((F_by_Fleet[2]/(F_by_Fleet[2]+F_by_Fleet[3]+F_by_Fleet[4]))/(TimeFit$`F:_1`[TimeFit$Yr==i]/(TimeFit$`F:_1`[TimeFit$Yr==i]+TimeFit$`F:_2`[TimeFit$Yr==i]+TimeFit$`F:_3`[TimeFit$Yr==i]))-1)+1
        if(abs(1-Fmult3F1)<=Relative.F.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        Fmult3F2<-rand2*((F_by_Fleet[3]/(F_by_Fleet[2]+F_by_Fleet[3]+F_by_Fleet[4]))/(TimeFit$`F:_2`[TimeFit$Yr==i]/(TimeFit$`F:_1`[TimeFit$Yr==i]+TimeFit$`F:_2`[TimeFit$Yr==i]+TimeFit$`F:_3`[TimeFit$Yr==i]))-1)+1
        if(abs(1-Fmult3F2)<=Relative.F.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        Fmult3F3<-rand2*((F_by_Fleet[4]/(F_by_Fleet[2]+F_by_Fleet[3]+F_by_Fleet[4]))/(TimeFit$`F:_3`[TimeFit$Yr==i]/(TimeFit$`F:_1`[TimeFit$Yr==i]+TimeFit$`F:_2`[TimeFit$Yr==i]+TimeFit$`F:_3`[TimeFit$Yr==i]))-1)+1
        if(abs(1-Fmult3F3)<=Relative.F.Threshold & keepFitting==FALSE){keepFitting<-FALSE}else{keepFitting<-TRUE}
        
      }
      #Finally update all the forecast F's by multiplying by all the scaling values calculated above.
      forecast$ForeCatch[(row+1),4]<-forecast$ForeCatch[(row+1),4]*min(1.5,max(0.5,DepletionScale*Fmult1*Fmult2G1*Fmult3F1))
      forecast$ForeCatch[(row+2),4]<-forecast$ForeCatch[(row+2),4]*min(1.5,max(0.5,DepletionScale*Fmult1*Fmult2G2*Fmult3F2))
      forecast$ForeCatch[(row+3),4]<-forecast$ForeCatch[(row+3),4]*min(1.5,max(0.5,DepletionScale*Fmult1*Fmult2G2*Fmult3F3))
      
      #Output scaling values to console so user can see if the model is converging
      print(paste(DepletionScale,Fmult1,Fmult2G1,Fmult2G2,Fmult3F1,Fmult3F2,Fmult3F3))
    }
  
    #Write the modified forecast data out to a file and rerun projections
    SS_writeforecast(mylist=forecast,overwrite = TRUE)
    shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
    #If all values have converged check if this is the OFL loop or the P* loop
    if(keepFitting==FALSE){
      if(fitting_depletion==TRUE){
        #If in the OLF loop then reset to keep fitting and change the target from OFL to P*
        keepFitting<-TRUE
        fitting_depletion<-FALSE
        #Write out the OFL results to a new folder (replace any old folder that exists)
        if(dir.exists(paste0(assessment_dir,"/Depletion_target"))){
          unlink(paste0(assessment_dir,"/Depletion_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/Depletion_target"))
        temp.files<-list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/Depletion_target/",temp.files))
      }else{
        #If on the P* loop write out the P* results to a new folder and end the procedure
        if(dir.exists(paste0(assessment_dir,"/Pstar_target"))){
          unlink(paste0(assessment_dir,"/Pstar_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/Pstar_target"))
        temp.files<-list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/Pstar_target/",temp.files))
      }
    }
  }
}