##---------------------------------------------------------------------##
#                                                                       #
#  SS projection package to achieve target equilibrium stock status,    #
#  constant annual fishing mortality rate, fixed group catch            #
#  allocations, and historic relative fishing effort between fleets     #
#  within groups. The package is able to calculate OFL, ABC represented #
#  as a fraction of the OFL F, and F rebuild represented as a fixed     #
#  reduction in F during the specified rebuild period if needed.        #
#                                                                       #
##---------------------------------------------------------------------##
#                                                                       #
#  Author: Nathan Vaughan                                               #
#  Last update date: 8/23/21                                             #
#  Contact: nathan.vaughan@noaa.gov                                     #
#                                                                       #
##---------------------------------------------------------------------##

# This is now a generalized function to run projections that should work for most species. 
# The function is now designed to mostly read out settings from the forecast file the 
# same way you would generally set things up for a base SS projection.
# The only additional inputs needed are a proportion for F.ABC relative to F.OFL if 
# ABC projections are desired and a rebuild target year if rebuilding projections are desired.
#

run.projections<-function(assessment_dir, #Here you set the location of a previously fit SS stock assessment to perform projections
                          ABC_Fraction = NULL, #Set the ABC target as a fraction of the OFL target if NULL will not fit ABC projections
                          Rebuild_yr = NULL, #Set the rebuild target year if NULL will not fit rebuild projections
                          F_max = FALSE, #If true and forecast method is Fmsy will replace Fmsy with Fmax search (maximum yield per recruit)
                          Depletion.Threshold = 0.001, # These are all just thresholds for when 
                          Annual.F.Threshold = 0.001, # targets are acceptably achieved these default to a .1% change
                          Allocation.Threshold = 0.001, # increase them if run is too slow or reduce to improve fit if run is fast.
                          Step.Threshold = 0.001,
                          rec_devs = rep(0,100), #Input for custom rec_devs and below implementation error needs a vector of 100 values one 
                          Fcast_impl_error = rep(0,100), #for each year of forecast. Defaults to no rec devs or implementation error.  
                          OFL_complete = FALSE #Set to true if you already have a fit OFL run but want to add or adjust an ABC or Rebuild run
                          ) 
{
  library(r4ss)
  #Set large max print to avoid issues with writing out a large forecast file.
  options(max.print = 1000000)
  setwd(paste0(assessment_dir))
  
  #Read in all the model files and results 
  start <- SS_readstarter()
  dat <- SS_readdat(file = start$datfile, version = 3.3)
  ctl <- SS_readctl(file = start$ctlfile, version = 3.3, use_datlist = TRUE, datlist = dat)
  results <- SS_output(dir = getwd(), covar = FALSE, checkcor = FALSE)
  forecast <- SS_readforecast() 
  parlist <- SS_readpar_3.30(parfile = "ss.par", datsource = dat, ctlsource = ctl)
  
  #First set up a working director for running projections in (to avoid overwriting the base files with a failed model run)
  #then copy all of the assessment files to this working folder (ignore any output directories that have been previously created)
  if(OFL_complete == FALSE){
  if(dir.exists(paste0(assessment_dir,"/Working_dir"))){
    unlink(paste0(assessment_dir,"/Working_dir"), recursive = TRUE)
  }
  dir.create(paste0(assessment_dir,"/Working_dir"))
  temp.files <- list.files(path=assessment_dir)
  folders <- c(which(temp.files=="OFL_target"), which(temp.files=="ABC_target"), which(temp.files=="Rebuild_target"), which(temp.files=="Working_dir"))
  if(length(folders)>0){
    temp.files <- temp.files[-folders]
  }
  file.copy(from = paste0(assessment_dir,'/',temp.files), to = paste0(assessment_dir,"/Working_dir/",temp.files))
  }else{
    temp.files <- list.files(path=paste0(assessment_dir,"/Working_dir"))
    folders <- c(which(temp.files=="OFL_target"), which(temp.files=="ABC_target"), which(temp.files=="Rebuild_target"), which(temp.files=="Working_dir"))
    if(length(folders)>0){
      temp.files <- temp.files[-folders]
    }
    unlink(paste0(assessment_dir,"/Working_dir/",temp.files))
    temp.files <- list.files(path=paste0(assessment_dir,"/OFL_target"))
    file.copy(from = paste0(assessment_dir,'/OFL_target/',temp.files), to = paste0(assessment_dir,"/Working_dir/",temp.files))
  }
  
  #Set the new working directory 
  setwd(paste0(assessment_dir,"/Working_dir"))
  
  #Modify assessment files to produce expected timeseries length and outputs
  #100 year projections to allow equilibrium to be achieved
  forecast$Nforecastyrs <- 100
  
  #Need to set recdevs and implementation error for all projection years 
  #so that reading from par file is possible
  parlist$recdev_forecast <- matrix(NA, nrow = 100, ncol = 2)
  parlist$Fcast_impl_error <- matrix(NA, nrow = 100, ncol = 2)
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
  
  #Identify the column numbers for Catch, F, SSB, Recruits, etc
  Catch_cols <- grep("retain(B)", names(TimeFit), fixed = TRUE)
  F_cols <- grep("F", names(TimeFit), fixed = TRUE)
  
  TimeFit2 <- aggregate(TimeFit[,sort(c(2,4,7,Catch_cols,F_cols))],by=list(TimeFit$Yr,TimeFit$Seas),FUN=sum)[,-c(3,4,5)]
  names(TimeFit2)[c(1,2)] <- c("Yr", "Seas")
  
  Catch_cols2 <- grep("retain(B)", names(TimeFit2), fixed = TRUE)
  F_cols2 <- grep("F", names(TimeFit2), fixed = TRUE)
  
  TimeFit3 <- aggregate(TimeFit[,sort(c(2,7,8,Catch_cols,F_cols))],by=list(TimeFit$Yr),FUN=sum)[,-2]
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
  
  #Build a projection forcast matrix of F values by fleet/season/year which will adjusted to achieve target stock status, F, and allocations.
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
  
  
  #Extract fixed forecast values from the forecast file these will be fixed at these values for the projections
  #This is used to implement recent catches or fixed harvest from an independent fleet such as shrimp bycatch.
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
  
  #Replace the temporary average F values with specified fixed inputs
  forecast_F[fixed_ref,c(4,5)]<-Fixed_catch_target[,c(4,5)]
  
  #This reference identifies inputs subject to a rebuilding period F 
  rebuild_ref <- which(forecast_F[,1]<=Rebuild_yr)
  
  #Here all input values are assigned to an allocation group if needed and
  #relative landings targets are identified
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
  
  keepFitting <- TRUE
  loop <- 0
  subloop <- 0
  
  F_adjust1 <- F_adjust2 <- 1
  F_adjust3 <- rep(1,100*length(seasons)*length(F_cols))
  search_step <- 0.1 
  Fmult1 <- Fmult2 <- Fmult3 <- rep(1.01,100*length(seasons)*length(F_cols))
  Fmult2a <- Fmult2b <- 1
  First_run<-TRUE
  
  if(OFL_complete == FALSE){
    #Save all the modified files and then perform a base run of SS so that output is specified correctly with 
    #a 100 year projection series.
    
    SS_writepar_3.30(parlist = parlist,outfile="ss.par",overwrite = TRUE)
    SS_writeforecast(mylist=forecast,overwrite = TRUE)
    SS_writestarter(mylist=start,overwrite = TRUE)
    
    shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
    
    #Begin the search in the OFL phase
    fitting_OFL <- TRUE
    fitting_ABC <- FALSE
    fitting_Rebuild <- FALSE
    
    MSY.Fit <- data.frame(catch=c(0),Ave.F=c(0),depletion=c(0),target.depletion=c(0))
    method <- "OFL"
    
  
  }else{
    #Save all the modified files and then set search to begin in ABC or Rebuild phase
    
    start <- SS_readstarter()
    dat <- SS_readdat(file = start$datfile, version = 3.3)
    ctl <- SS_readctl(file = start$ctlfile, version = 3.3, use_datlist = TRUE, datlist = dat)
    results <- SS_output(dir = getwd(), covar = FALSE, checkcor = FALSE)
    forecast <- SS_readforecast() 
    parlist <- SS_readpar_3.30(parfile = "ss.par", datsource = dat, ctlsource = ctl)
    
    fitting_OFL <- FALSE
    if(!is.null(ABC_Fraction)){
      fitting_ABC <- TRUE
      fitting_Rebuild <- FALSE
      method <- "ABC"
    }else if(!is.null(Rebuild_yr)){
      fitting_ABC <- FALSE
      fitting_Rebuild <- TRUE
      method <- "Rebuild"
    }else{
      keepFitting <- FALSE
      fitting_ABC <- FALSE
      fitting_Rebuild <- FALSE
      method <- NULL
    }
  }
  
  #Set up plot window for production of search diagnostic plots depending on target specifications
  #These plots were largely for diagnostic testing during code development but can allow you to see what
  #is going wrong if a future bug does occur.
  if(Forecast_target==2 & fitting_OFL==TRUE){
    par(mfrow=c(5,2))
  }else{
    par(mfrow=c(4,2))
  }
  
  #Now start a loop of projecting and modifying fixed F's until the desired 
  #landings projections are achieved
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
    
    TimeFit3 <- aggregate(TimeFit[,sort(c(2,7,8,Catch_cols,F_cols))],by=list(TimeFit$Yr),FUN=sum)[,-2]
    names(TimeFit3)[c(1)] <- c("Yr")
    
    #If reading in results of a previous OFL run then set the F_OFL and F.ABC values before begining ABC/Rebuild loops
    if(OFL_complete==TRUE & First_run==TRUE){
      F_report<-SPRfit$F_report
      FScale<-median(F_report[(length(F_report)-89):length(F_report)])
      F_OFL<-FScale
      if(!is.null(ABC_Fraction)){
        F.ABC<-ABC_Fraction*FScale
      }else{
        F.ABC<-FScale
      }
      First_run<-FALSE
    }
    
    if(fitting_OFL==TRUE){
      
      #Calculate the average F at equilibrium that all F's will be scaled to in order
      #to achieve equal F in every year. As depletion approaches the target value this 
      #F will approach F(OFL).
      F_report<-SPRfit$F_report
      FScale<-median(F_report[(length(F_report)-89):length(F_report)])
      F_OFL<-FScale
      if(!is.null(ABC_Fraction)){
        F.ABC<-ABC_Fraction*FScale
      }else{
        F.ABC<-FScale
      }
      
      #Calculate depletion target adjustment scale depending on the specified target (SPR ratio, SSB ratio, or true MSY)   
      if(Forecast_target==1){
        search_step<-0.00001
        Target.Depletion <- forecast[["SPRtarget"]]
        Depletion<-SPRfit$SPR
        
        Achieved.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        
        DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
        
        Depletion_R<-TimeFit3$SpawnBio/Virgin_bio
        Target.Rebuild <- mean(Depletion_R[(length(Depletion_R)-9):length(Depletion_R)])
        
      }else if(Forecast_target==2){
        Depletion <- TimeFit3$SpawnBio/Virgin_bio
        Achieved.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
        if(First_run == TRUE){
          Target.Depletion <- Achieved.Depletion
          First_run <- FALSE
        }
        Target.Rebuild <- Target.Depletion
        
        if(max(abs(1-Fmult3))>Allocation.Threshold | max(abs(1-Fmult2))>Annual.F.Threshold | max(abs(1-Fmult1))>Depletion.Threshold){
          loop<-loop-1
          subloop<-subloop+1
          if(F_max==TRUE){
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),"Recruit_0"])
          }else{
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])
          }
          MSY.Fit[1,] <- c(Achieved.Catch,FScale,Achieved.Depletion,Target.Depletion)
        }else{
          subloop<-0
          if(F_max==TRUE){
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),"Recruit_0"])
          }else{
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])
          }
          MSY.Fit <- rbind(MSY.Fit[1,],MSY.Fit)
          MSY.Fit[1,] <- c(Achieved.Catch,FScale,Achieved.Depletion,Target.Depletion)
          if(loop>1){
            if(Achieved.Catch<Last_Achieved_Catch){
              search_step <- -0.5*search_step
            }
          
            Target.Depletion <- Target.Depletion+search_step
            
            min_diff <- which(abs(MSY.Fit[,4]-Target.Depletion)<0.001)
            if(length(min_diff)>0){
              Old.Catch <- MSY.Fit[min_diff[1],1]
              if(Old.Catch<Achieved.Catch){
                search_step <- -0.5*search_step
              }
              Target.Depletion <- Target.Depletion+search_step
              Achieved.Catch <- Old.Catch
            }
          }else{
            steps <- seq(0.1,0.9,0.1)
            New.Target.Depletion <- steps[which(abs(steps-Target.Depletion)==min(abs(steps-Target.Depletion)))[1]]
            if(New.Target.Depletion<Target.Depletion){
              search_step <- -1*search_step
            }
            Target.Depletion <- New.Target.Depletion
          }
          Last_Achieved_Catch <- Achieved.Catch
        }
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        if(F_max==TRUE){
          plot(x=TimeFit3[,"Yr"],y=apply(TimeFit3[,Catch_cols3],1,sum)/TimeFit3[,"Recruit_0"],xlab="year",ylab="Total Yield Per Recruit",main = paste0(method," loop = ",loop,".",subloop))
          plot(x=MSY.Fit[,"depletion"],y=MSY.Fit[,"catch"],xlim=c(0.9*min(MSY.Fit[,3:4]),1.1*max(MSY.Fit[,3:4])),xlab="Esimated depletion",ylab="Total Yield Per Recruit",main = paste0(method," loop = ",loop,".",subloop))
          lines(x=c(MSY.Fit[1,c(4,4)]),y=c(0,2*max(MSY.Fit[,"catch"])),col="dark red")
          points(x=MSY.Fit[1,3],y=MSY.Fit[1,1],pch=16,col="dark blue")
        }else{
          plot(x=TimeFit3[,"Yr"],y=apply(TimeFit3[,Catch_cols3],1,sum),xlab="year",ylab="Total Yield",main = paste0(method," loop = ",loop,".",subloop))
          plot(x=MSY.Fit[,"depletion"],y=MSY.Fit[,"catch"],xlim=c(0.9*min(MSY.Fit[,3:4]),1.1*max(MSY.Fit[,3:4])),xlab="Esimated depletion",ylab="Total Yield",main = paste0(method," loop = ",loop,".",subloop))
          lines(x=c(MSY.Fit[1,c(4,4)]),y=c(0,2*max(MSY.Fit[,"catch"])),col="dark red")
          points(x=MSY.Fit[1,3],y=MSY.Fit[1,1],pch=16,col="dark blue")
        }
      }else if(Forecast_target==3){
        search_step<-0.00001
        Target.Depletion <- Target.Rebuild <- forecast[["Btarget"]]
        Depletion<-TimeFit3$SpawnBio/Virgin_bio

        Achieved.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
      }
    }else if(fitting_ABC==TRUE){
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the OFL MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that OFL search is complete
      FScale<-F.ABC #Set the F target to F.ABC for rescaling annual F values
    }else if(fitting_Rebuild==TRUE){
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the OFL MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that OFL search is complete
      FScale<-F_OFL #Set the F target to F_OFL for rescaling annual F values in years after the rebuild period.
      F_report<-SPRfit$F_report
      F_Rebuild_Scale<-F_report[SPRfit$Yr==Rebuild_yr]
      Depletion<-TimeFit3$SpawnBio/Virgin_bio
      Achieved.Rebuild <- mean(Depletion[SPRfit$Yr==Rebuild_yr])
      Rebuild.Scale <- (1-Target.Rebuild)/(1-Achieved.Rebuild)
      Rebuild.Scale <- (-log(1-((1-exp(-F_Rebuild_Scale))*Rebuild.Scale))/F_Rebuild_Scale)
      Rebuild.Scale <- Rebuild.Scale*F_Rebuild_Scale
      Rebuild.Scale <- min(Rebuild.Scale,FScale)
    }
    
    #Fmult2 calculations define the multiplier for adjusting annual F values
    #Zero catch years are identified first to prevent divide by zero errors in the scaling and
    #allow tell the search the target has been achieved
    zero_catch <- which(SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))]==0)
    if(length(zero_catch)>0){
      Fmult2[zero_catch] <- 1
      Fmult2[-zero_catch] <- FScale/SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))][-zero_catch]
    }else{
      Fmult2 <- FScale/SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))]
    }
    #If in a rebuild search phase the rebuild years are now adjusted independently of the later F_OFL years
    if(fitting_Rebuild==TRUE){
      Fmult2[rebuild_ref] <- Rebuild.Scale/SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))][rebuild_ref]
    }
    
    #Here 
    if((loop>1 |  subloop>2) & DepletionScale>0 & fitting_OFL==TRUE & (Forecast_target!=2 | subloop>2)){
      F_adjust1 <- (F_adjust1 + 1)/2
      F_adjust1 <- F_adjust1*(Last_Mult1-1)/(Last_Mult1-DepletionScale)
      
      if(F_adjust1>0){
        if(DepletionScale<1){
          DepletionScale <- exp(log(DepletionScale)*F_adjust1)
        }
        if(DepletionScale>1){
          DepletionScale <- ((DepletionScale-1)*F_adjust1+1)
        }
      }else{
        F_adjust1 <- 1
      }
    }else if(loop>1 & min(Fmult2)>0 & fitting_ABC==TRUE){
      F_adjust2 <- (F_adjust2 + 1)/2
      F_adjust2 <- F_adjust2*(Last_Mult2-1)/(Last_Mult2-mean(Fmult2[-fixed_ref]))
      if(F_adjust2 > 0){
        if(length(Fmult2[Fmult2<1])>0){
          Fmult2[Fmult2<1] <- exp(log(Fmult2[Fmult2<1])*F_adjust2)
        }
        if(length(Fmult2[Fmult2>1])>0){
          Fmult2[Fmult2>1] <- ((Fmult2[Fmult2>1]-1)*F_adjust2+1)
        }
      }else{
        F_adjust2 <- 1
      }
    }else if(loop>1 & min(Fmult2)>0 & fitting_Rebuild==TRUE){
      F_adjust2a <- (Last_Mult2a-1)/(Last_Mult2a-mean(Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]]))
      F_adjust2b <- (Last_Mult2b-1)/(Last_Mult2b-mean(Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))]))
      if(!is.na(F_adjust2a)){
      if(F_adjust2a > 0){
        if(exists("Fmult2a")){
          Fmult2a <- (Fmult2a + 1)/2
        }
        Fmult2a <- Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]]
        if(length(Fmult2a[Fmult2a<1])>0){
          Fmult2a[Fmult2a<1] <- exp(log(Fmult2a[Fmult2a<1])*F_adjust2a)
        }
        if(length(Fmult2a[Fmult2a>1])>0){
          Fmult2a[Fmult2a>1] <- ((Fmult2a[Fmult2a>1]-1)*F_adjust2a+1)
        }
        Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]] <- Fmult2a
      }else{
        F_adjust2a <- 1
      }}else{
        F_adjust2a <- 1
      }
      if(!is.na(F_adjust2b)){
      if(F_adjust2b > 0){
        if(exists("Fmult2b")){
          Fmult2b <- (Fmult2b + 1)/2
        }
        Fmult2b <- Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))]
        if(length(Fmult2b[Fmult2b<1])>0){
          Fmult2b[Fmult2b<1] <- exp(log(Fmult2b[Fmult2b<1])*F_adjust2b)
        }
        if(length(Fmult2b[Fmult2b>1])>0){
          Fmult2b[Fmult2b>1] <- ((Fmult2b[Fmult2b>1]-1)*F_adjust2b+1)
        }
        Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))] <- Fmult2b
      }else{
        F_adjust2b <- 1
      }}else{
        F_adjust2b <- 1
      }
    }
    
    Fmult1 <- rep(DepletionScale,100*length(seasons)*length(F_cols))
    
    if(FScale > 0){               
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
      Fmult3 <- (0.5*(Allocations[,5]/Allocations[,6]-1)+1)
    }else{
      Fmult3 <- rep(1,100*length(seasons)*length(F_cols))
    }
    Fmult1[fixed_ref] <- 1
    Fmult2[fixed_ref] <- 1
    Fmult3[fixed_ref] <- 1
    Comb_Mult <- Fmult1*Fmult2*Fmult3
    #Comb_Mult <- pmin(pmax(0.2,Comb_Mult),2)
                   
    Last_Mult1 <- DepletionScale
    Last_Mult2 <- mean(Fmult2[-fixed_ref])
    Last_Mult2a <- mean(Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]])
    Last_Mult2b <- mean(Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))])
    
    col_options <- c("black","dark red","dark green","dark blue","orange","purple","red","green","blue","brown","pink","yellow",colors())
    point_options <- c(16,15,17,18,8,9,10,11,12,13,0,1,2,3,4,5,6,14,21,22,23,24,25,19,20)
    plot(Fmult1,xlab="year/season/fleet",ylab="Depletion Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(rep(F_adjust1,100*length(seasons)*length(F_cols)),xlab="year/season/fleet",ylab="Depletion Optimization Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(Fmult2,xlab="year/season/fleet",ylab="F Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(rep(F_adjust2,100*length(seasons)*length(F_cols)),xlab="year/season/fleet",ylab="F Optimization Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(Fmult3,xlab="year/season/fleet",ylab="Allocation Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    plot(F_adjust3,xlab="year/season/fleet",ylab="Allocation Optimization Adjustment",col=rep(col_options[seq_along(F_cols)],100*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),100),main = paste0(method," loop = ",loop))
    
    if(max(abs(1-Fmult1))>=Depletion.Threshold | max(abs(1-Fmult2))>=Annual.F.Threshold | max(abs(1-Fmult3))>=Allocation.Threshold | abs(search_step)>Step.Threshold){keepFitting<-TRUE}else{keepFitting<-FALSE}
    if(FScale==0 & loop>2){keepFitting<-FALSE}
    forecast_F[,4] <- forecast_F[,4]*Comb_Mult
    forecast_F[fixed_ref,4] <- Fixed_catch_target[,4]
    forecast[["ForeCatch"]] <- forecast_F
    #Write the modified forecast data out to a file and rerun projections
    unlink(paste0(getwd(),"/forecast.ss"))
    SS_writeforecast(mylist=forecast,overwrite = TRUE)
    shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
    #If all values have converged check if this is the OFL loop or the P* loop
    if(keepFitting==FALSE){
      if(fitting_OFL==TRUE){
        fitting_OFL <- FALSE
        par(mfrow=c(4,2))
        if(!is.null(ABC_Fraction)){
          #If in the OLF loop then reset to keep fitting and change the target from OFL to P*
          loop <- 0
          method <- "ABC"
          keepFitting <- TRUE
          fitting_ABC <- TRUE
        }else if(!is.null(Rebuild_yr)){
          #If in the OLF loop then reset to keep fitting and change the target from OFL to Rebuild
          loop <- 0
          method <- "Rebuild"
          keepFitting <- TRUE
          fitting_Rebuild <- TRUE
        }
        #Write out the OFL results to a new folder (replace any old folder that exists)
        if(dir.exists(paste0(assessment_dir,"/OFL_target"))){
          unlink(paste0(assessment_dir,"/OFL_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/OFL_target"))
        temp.files <- list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/OFL_target/",temp.files))
      }else if(fitting_ABC==TRUE){
        fitting_ABC <- FALSE
        if(!is.null(Rebuild_yr)){
          #If in the ABC loop then reset to keep fitting and change the target from P* to Rebuild
          loop <- 0
          method <- "Rebuild"
          keepFitting <- TRUE
          fitting_Rebuild <- TRUE
        }
        #If on the ABC loop write out the P* results to a new folder (replace any old folder that exists)
        if(dir.exists(paste0(assessment_dir,"/ABC_target"))){
          unlink(paste0(assessment_dir,"/ABC_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/ABC_target"))
        temp.files <- list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/ABC_target/",temp.files))
      }else if(fitting_Rebuild==TRUE){
        fitting_Rebuild <- FALSE
        #If on the Rebuild loop write out the Rebuild results to a new folder (replace any old folder that exists) and end the fitting process
        if(dir.exists(paste0(assessment_dir,"/Rebuild_target"))){
          unlink(paste0(assessment_dir,"/Rebuild_target"),recursive = TRUE)
        }
        dir.create(paste0(assessment_dir,"/Rebuild_target"))
        temp.files <- list.files(path=paste0(assessment_dir,"/Working_dir"))
        file.copy(from=paste0(assessment_dir,'/Working_dir/',temp.files),to=paste0(assessment_dir,"/Rebuild_target/",temp.files))
      }
    }
  }
}