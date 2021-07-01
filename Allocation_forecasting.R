remotes::install_github(repo="r4ss/r4ss",ref="development")

library(r4ss)

#Test 

#This is currently set up for yellowtail snapper. Will need work to generalize to all assessment formats.

run.projections(assessment_dir="C:/Users/Nathan/Desktop/run25/run25.10")

run.projections<-function(assessment_dir, #Here you set the location of a previously fit SS stock assessment
                          PstarFraction = 0.75, #Set the P* target as a fraction of the OFL target
                          Depletion.Threshold = 0.0001, # These are all just thresholds for when 
                          Annual.F.Threshold = 0.0001, # targets are acceptably achieved these default to a .01% change
                          Allocation.Threshold = 0.0001, # increase them if run is too slow.
                          rec_devs = rep(0,100), 
                          Fcast_impl_error = rep(0,100) #defaults to no rec devs or implementation error  
                          ) 
{
  #Set large max print to avoid issues with writing out a large forecast file.
  options(max.print = 1000000)
  setwd(paste0(assessment_dir))
  #First set up a working director for running projections in (to avoid overwriting the base files with a failed model run)
  #then copy all of the assessment files to this working folder (ignore any output directories that have been previously created)
  if(dir.exists(paste0(assessment_dir,"/Working_dir"))){
    unlink(paste0(assessment_dir,"/Working_dir"), recursive = TRUE)
  }
  dir.create(paste0(assessment_dir,"/Working_dir"))
  temp.files <- list.files(path=assessment_dir)
  folders <- c(which(temp.files=="Depletion_target"), which(temp.files=="Pstar_target"), which(temp.files=="Working_dir"))
  if(length(folders)>0){
    temp.files <- temp.files[-folders]
  }
  file.copy(from = paste0(assessment_dir,'/',temp.files), to = paste0(assessment_dir,"/Working_dir/",temp.files))
  
  #Set working directory and read in the assessment files
  setwd(paste0(assessment_dir,"/Working_dir"))
  start <- SS_readstarter()
  dat <- SS_readdat(file = start$datfile, version = 3.3)
  ctl <- SS_readctl(file = start$ctlfile, version = 3.3, use_datlist = TRUE, datlist = dat)
  results <- SS_output(dir = getwd(), covar = FALSE, checkcor = FALSE)
  forecast <- SS_readforecast() 
  parlist <- SS_readpar_3.30(parfile = "ss.par", datsource = dat, ctlsource = ctl)
  
  
  
  #Modify assessment files to produce expected timeseries length and outputs
  #100 year projections allow equilibrium to be achieved
  forecast$Nforecastyrs <- 100
  parlist$recdev_forecast <- matrix(NA, nrow = 100, ncol = 2)
  parlist$Fcast_impl_error <- matrix(NA, nrow = 100, ncol = 2)
  #Need to set recdevs and implementation error for all projection years 
  #so that reading from par file is possible
  parlist$recdev_forecast[,1] <- (dat[["endyr"]]+1):(dat[["endyr"]]+100)
  parlist$Fcast_impl_error[,1] <- (dat[["endyr"]]+1):(dat[["endyr"]]+100)
  parlist$recdev_forecast[,2] <- rec_devs
  parlist$Fcast_impl_error[,2] <- Fcast_impl_error
  colnames(parlist$recdev_forecast) <- c("year","recdev")
  colnames(parlist$Fcast_impl_error) <- c("year","impl_error")
  
  #Adjust the starter file to read from par file, perform no fitting (This should already have been done),
  #and set the depletion value to be relative to unexploited biomass and have no scaling 
  #(so that correct depletion target can be found).
  start$init_values_src <- 1
  start$last_estimation_phase <- 0
  start$depl_basis <- 1
  start$depl_denom_frac <- 1
  start$SPR_basis <- 4
  start$F_report_units <- 1
  start$F_report_basis <- 0
  
  
  #Get the timeseries of historic/projected catches and F
  TimeFit <- results$timeseries
  #TimeFit <- TimeFit[TimeFit[,"Yr"]>dat[["endyr"]],]
  
  Catch_cols <- grep("retain(B)", names(TimeFit), fixed = TRUE)
  F_cols <- grep("F", names(TimeFit), fixed = TRUE)
  
  TimeFit2 <- aggregate(TimeFit[,sort(c(2,4,7,Catch_cols,F_cols))],by=list(TimeFit$Yr,TimeFit$Seas),FUN=sum)[,-c(3,4,5)]
  names(TimeFit2)[c(1,2)] <- c("Yr", "Seas")
  
  Catch_cols2 <- grep("retain(B)", names(TimeFit2), fixed = TRUE)
  F_cols2 <- grep("F", names(TimeFit2), fixed = TRUE)
  
  TimeFit3 <- aggregate(TimeFit[,sort(c(2,7,Catch_cols,F_cols))],by=list(TimeFit$Yr),FUN=sum)[,-2]
  names(TimeFit3)[c(1)] <- c("Yr")
  Virgin_bio <- TimeFit3$SpawnBio[1]
  
  Catch_cols3 <- grep("retain(B)", names(TimeFit3), fixed = TRUE)
  F_cols3 <- grep("F", names(TimeFit3), fixed = TRUE)
  
  #Set all future projections to fish at constant apical F that matches recent years
  #get the years of timeseries F's based on the forecast year range
  TargetYears <- TimeFit2[TimeFit2$Yr>=forecast$Fcast_years[3] & TimeFit2$Yr<=forecast$Fcast_years[4],]
  TargetYears <- TargetYears[,c(2,F_cols2)]
  seasons <- unique(TargetYears[,1])
  F_by_Fleet_seas <- as.data.frame(matrix(apply(TargetYears[TargetYears[,1]==seasons[1],], 2, mean),nrow=1,ncol=(length(F_cols)+1)))
  if(length(seasons)>1){
    for(i in seasons[-1]){
      F_by_Fleet_seas <- rbind(F_by_Fleet_seas,apply(TargetYears[TargetYears[,1]==i,], 2, mean))
    }
  }
  Forecast_target <- forecast[["Forecast"]]
  if(!is.element(Forecast_target,c(1,2,3))){
    stop("forecast should be set to either 1, 2, or 3 so we know what the target is")
  }
  
  forecast_F<-matrix(1,nrow=100*length(seasons)*length(F_cols),ncol=5)
  forecast_F[,1]<-sort(rep((dat[["endyr"]]+1):(dat[["endyr"]]+100),length(seasons)*length(F_cols)))
  forecast_F[,2]<-rep(seasons,100*length(F_cols))
  forecast_F[,3]<-rep(1:length(F_cols),100*(length(seasons)))
  for(i in seasons){
    forecast_F[forecast_F[,2]==i,4]<-unlist(rep(F_by_Fleet_seas[F_by_Fleet_seas[,1]==i,-1],100))
    forecast_F[forecast_F[,2]==i,5]<-rep(99,length(F_cols)*100)
  }
  forecast_F<-as.data.frame(forecast_F)
  colnames(forecast_F)<-c("Year","Seas","Fleet","Catch or F","Basis")
  
  #Extract fixed forecast values from the forecast file these will be fixed at these values
  if(!is.null(forecast[["ForeCatch"]])){
    Fixed_catch_basis <- forecast[["InputBasis"]]
    Fixed_catch_target <- forecast[["ForeCatch"]]
    fixed_ref <- seq_along(Fixed_catch_target[,1])
    if(length(names(Fixed_catch_target))==4){
      new_names <- c(names(Fixed_catch_target),"Basis")
      Fixed_catch_target <- cbind(Fixed_catch_target,rep(Fixed_catch_basis,length(Fixed_catch_target[,1])))
      names(Fixed_catch_target) <- new_names
    }
    for(i in seq_along(Fixed_catch_target[,1])){
      fixed_ref[i]<-which(forecast_F[,1]==Fixed_catch_target[i,1] & 
                          forecast_F[,2]==Fixed_catch_target[i,2] & 
                          forecast_F[,3]==Fixed_catch_target[i,3])
    }
  }
  
  forecast_F[fixed_ref,c(4,5)]<-Fixed_catch_target[,c(4,5)]
  
  n_groups <- forecast[["N_allocation_groups"]]
  groups <- rep(0,length(F_cols))
  Allocations <- forecast_F[,c(1,2,3,4,5,5)]
  Allocations[,c(4)] <- 0
  Allocations[,c(5,6)] <- 1
  names(Allocations) <- c("Year","Seas","Fleet","Group","Target","Achieved")
  if(n_groups>0){
    for(i in seq_along(forecast[["fleet_assignment_to_allocation_group"]][,"Fleet"])){
      groups[forecast[["fleet_assignment_to_allocation_group"]][i,"Fleet"]] <- forecast[["fleet_assignment_to_allocation_group"]][i,"Group"]
      Allocations[Allocations[,"Fleet"]==forecast[["fleet_assignment_to_allocation_group"]][i,"Fleet"],4] <- forecast[["fleet_assignment_to_allocation_group"]][i,"Group"]
    }
    
    alloc <- forecast[["allocation_among_groups"]][order(forecast[["allocation_among_groups"]][,"Year"]),]
    for(i in seq_along(alloc[,1])){
      for(j in seq_along(alloc[i,-1])){
        Allocations[Allocations[,1]>=alloc[,"Year"] & Allocations[,4]==j,5] <- alloc[i,(j+1)]/sum(alloc[i,-1])
      }
    }
  }
  
  forecast[["Forecast"]] <- 4
  forecast[["InputBasis"]] <- -1
  forecast[["ForeCatch"]] <- forecast_F
  forecast[["FirstYear_for_caps_and_allocations"]] <- (dat[["endyr"]]+101)
  #Save all the modified files and then perform a base run of SS so that output is dimensioned correctly with 
  #a 100 year projection series
  SS_writepar_3.30(parlist = parlist,outfile="ss.par",overwrite = TRUE)
  SS_writeforecast(mylist=forecast,overwrite = TRUE)
  SS_writestarter(mylist=start,overwrite = TRUE)
  
  shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
  
  #Now start a loop of projecting and modifying fixed F's until the desired 
  #catch projections are achieved
  keepFitting <- TRUE
  fitting_depletion <- TRUE
  loop <- 0
  MSY.Fit <- data.frame(catch=c(0),Ave.F=c(0),weight=c(1))
  method <- "MSY/proxy target"
  par(mfrow=c(3,2))
  FAdjust <- matrix(c(rep(1,100*length(seasons)*length(F_cols)),
                      rep(0.1,100*length(seasons)*length(F_cols)),
                      rep(10,100*length(seasons)*length(F_cols))),
                    nrow=100*length(seasons)*length(F_cols),
                    ncol=3, byrow = FALSE)
  
  while(keepFitting){
    #Read in the SS results for landings and stock status to determine if desired
    #targets have been achieved
    resultsFit <- SS_output(dir=getwd(),covar=FALSE,checkcor = FALSE)
    TimeFit <- resultsFit[["timeseries"]]
    TimeFit <- TimeFit[TimeFit[,"Yr"]>dat[["endyr"]],]
    SPRfit <- resultsFit[["sprseries"]]
    SPRfit <- SPRfit[SPRfit[,"Yr"]>dat[["endyr"]],]
    loop <- loop + 1
    plot(SPRfit[SPRfit[,"Yr"]>=dat[["endyr"]],"F_report"],xlab="year",ylab="F",main = paste0(method," loop = ",loop))
    plot(SPRfit[SPRfit[,"Yr"]>=dat[["endyr"]],"Deplete"],xlab="year",ylab="Depletion",main = paste0(method," loop = ",loop))
    
    
    TimeFit2 <- aggregate(TimeFit[,sort(c(2,4,7,Catch_cols,F_cols))],by=list(TimeFit$Yr,TimeFit$Seas),FUN=sum)[,-c(3,4,5)]
    names(TimeFit2)[c(1,2)] <- c("Yr", "Seas")
    
    TimeFit3 <- aggregate(TimeFit[,sort(c(2,7,Catch_cols,F_cols))],by=list(TimeFit$Yr),FUN=sum)[,-2]
    names(TimeFit3)[c(1)] <- c("Yr")
    
    
    if(fitting_depletion==TRUE){
      
      #Calculate the average F at equilibrium that all F's will be scaled to in order
      #to achieve equal F in every year. As depletion approaches the target value this 
      #F will approach F(OFL).
      F_report<-SPRfit$F_report
      FScale<-mean(F_report[(length(F_report)-9):length(F_report)])
      Fpstar<-FScale
      
      if(Forecast_target==1){
        Target.Depletion <- forecast[["SPRtarget"]]
        Depletion<-SPRfit$SPR
        
        Achieved.Depletion <- mean(Depletion[(length(Depletion)-9):length(Depletion)])
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        #DepletionScale <- Achieved.Depletion/Target.Depletion
        
        DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
      }else if(Forecast_target==2){
        Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])
        MSY.Fit <- rbind(MSY.Fit[1,],MSY.Fit)
        MSY.Fit[1,] <- c(Achieved.Catch,FScale,1)
        MSY.Fit[,"weight"] <- MSY.Fit[,"catch"]/max(MSY.Fit[,"catch"])
        parab<-function(b,MSY.Fit){
          exp<-b[1]-b[2]*(b[3]-MSY.Fit[,2])^2
          resid<-sum(MSY.Fit[,3]*(MSY.Fit[,1]-exp)^2)
        }
        if(loop<10){
          f.msy <- runif(1,0.75,1.25)*MSY.Fit[MSY.Fit[,"catch"]==max(MSY.Fit[,"catch"]),"Ave.F"][1]
        }else{
          if(loop==10){
            b <- c(max(MSY.Fit[,"catch"]),(max(MSY.Fit[,"catch"])/((MSY.Fit[MSY.Fit[,"catch"]==max(MSY.Fit[,"catch"]),"Ave.F"][1])^2)),MSY.Fit[MSY.Fit[,"catch"]==max(MSY.Fit[,"catch"]),"Ave.F"][1])
          }
          msy.estim <- optim(par=b,fn=parab,MSY.Fit=MSY.Fit)
          b <- msy.estim$par
          c.msy <-  msy.estim$par[1]
          f.msy <-  msy.estim$par[3]
        }
        DepletionScale <- f.msy/FScale
      }else if(Forecast_target==3){
        Target.Depletion <- forecast[["Btarget"]]
        Depletion<-TimeFit3$SpawnBio/Virgin_bio

        Achieved.Depletion <- mean(Depletion[(length(Depletion)-9):length(Depletion)])
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        #DepletionScale <- Achieved.Depletion/Target.Depletion
        
        DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
      }
    }else{
      DepletionScale<-1
      FScale<-PstarFraction*Fpstar
    }
    
    if(n_groups>0){
      Catch_temp <- TimeFit3[,Catch_cols3]
      Catch_tot <- apply(Catch_temp[,which(groups!=0)],1,sum)
      for(i in 1:n_groups){
        sort.mat <- matrix(NA, nrow = 100*length(seasons)*length(which(groups==i)), ncol = 2)
        sort.mat[,1] <- rep(1:100,length(seasons)*length(which(groups==i)))
        sort.mat[,2] <- rep(apply(Catch_temp[,which(groups==i)],1,sum)/Catch_tot,length(seasons)*length(which(groups==i)))
        sort.mat <- sort.mat[order(sort.mat[,1]),]
        Allocations[Allocations[,4]==i,6] <- sort.mat[,2]
      }
    }
    
    Fmult1 <- rep(DepletionScale,100*length(seasons)*length(F_cols))
    Fmult2 <- FScale/SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))]
    Fmult3 <- (0.5*(Allocations[,5]/Allocations[,6]-1)+1)
    #Fmult3 <- (-log(1-((1-exp(-forecast_F[,4]))*Fmult3))/forecast_F[,4])
    Fmult1[fixed_ref] <- 1
    Fmult2[fixed_ref] <- 1
    Fmult3[fixed_ref] <- 1
    Comb_Mult <- Fmult1*Fmult2*Fmult3
    Final_Mult <- Comb_Mult
    
    col_options <- c("black","dark red","dark green","dark blue","orange","purple","red","green","blue","brown","pink","yellow",colors())
    point_options <- c(16,15,17,18,8,9,10,11,12,13,0,1,2,3,4,5,6,14,21,22,23,24,25,19,20)
    plot(Fmult1,xlab="year/season/fleet",ylab="Depletion Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(Fmult2,xlab="year/season/fleet",ylab="F Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(Fmult3,xlab="year/season/fleet",ylab="Allocation Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(FAdjust[,1],xlab="year/season/fleet",ylab="Optimization Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    
    if(max(abs(1-Fmult1))>=Depletion.Threshold | max(abs(1-Fmult2))>=Annual.F.Threshold | max(abs(1-Fmult3))>=Allocation.Threshold  | loop<10){keepFitting<-TRUE}else{keepFitting<-FALSE}
    
    forecast_F[,4] <- forecast_F[,4]*Final_Mult
    forecast_F[fixed_ref,4] <- Fixed_catch_target[,4]
    forecast[["ForeCatch"]] <- forecast_F
    #Write the modified forecast data out to a file and rerun projections
    SS_writeforecast(mylist=forecast,overwrite = TRUE)
    shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
    #If all values have converged check if this is the OFL loop or the P* loop
    if(keepFitting==FALSE){
      if(fitting_depletion==TRUE){
        method <- "P*"
        #If in the OLF loop then reset to keep fitting and change the target from OFL to P*
        keepFitting <- TRUE
        fitting_depletion <- FALSE
        #Write out the OFL results to a new folder (replace any old folder that exists)
        if(dir.exists(paste0(assessment_dir,"/Depletion_target"))){
          unlink(paste0(assessment_dir,"/Depletion_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/Depletion_target"))
        temp.files <- list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/Depletion_target/",temp.files))
      }else{
        #If on the P* loop write out the P* results to a new folder and end the procedure
        if(dir.exists(paste0(assessment_dir,"/Pstar_target"))){
          unlink(paste0(assessment_dir,"/Pstar_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/Pstar_target"))
        temp.files <- list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/Pstar_target/",temp.files))
      }
    }
  }
}