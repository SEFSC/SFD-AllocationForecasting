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
#  Last update date: 7/9/21                                             #
#  Contact: nathan.vaughan@noaa.gov                                     #
#                                                                       #
##---------------------------------------------------------------------##

# remotes::install_github(repo="r4ss/r4ss",ref="development")
# library(r4ss)
# run.projections(assessment_dir="C:/Users/Nathan/Desktop/run25/run25.10", ABC_Fraction=.75, Rebuild_yr=2029)


#
# This is now a generalized function to run projections that should work for most species. 
# The function is now designed to mostly read out settings from the forecast file the 
# same way you would generally set things up for a base SS projection.
# The only additional inputs needed are a proportion for F.ABC relative to F.OFL if 
# ABC projections are desired and a rebuild target year if rebuilding projections are desired.
#

run.projections<-function(assessment_dir, #Here you set the location of a previously fit SS stock assessment to perform projections
                          ABC_Fraction = NULL, #Set the ABC target as a fraction of the OFL target if NULL will not fit ABC projections
                          Rebuild_yr = NULL, #Set the rebuild target year if NULL will not fit rebuild projections
                          Depletion.Threshold = 0.0001, # These are all just thresholds for when 
                          Annual.F.Threshold = 0.0001, # targets are acceptably achieved these default to a .01% change
                          Allocation.Threshold = 0.0001, # increase them if run is too slow.
                          rec_devs = rep(0,100), #Input for custom rec_devs and below implementation error needs a vector of 100 values one 
                          Fcast_impl_error = rep(0,100) #for each year of forecast. Defaults to no rec devs or implementation error.  
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
  folders <- c(which(temp.files=="OFL_target"), which(temp.files=="ABC_target"), which(temp.files=="Rebuild_target"), which(temp.files=="Working_dir"))
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
  
  rebuild_ref <- which(forecast_F[,1]<=Rebuild_yr)
  
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
  fitting_OFL <- TRUE
  fitting_ABC <- FALSE
  fitting_Rebuild <- FALSE
  loop <- 0
  MSY.Fit <- data.frame(catch=c(0),Ave.F=c(0),weight=c(1))
  method <- "MSY/proxy target"
  par(mfrow=c(4,2))
  F_adjust1 <- F_adjust2 <- 1
  F_adjust3 <- rep(1,100*length(seasons)*length(F_cols))
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
    
    
    if(fitting_OFL==TRUE){
      
      #Calculate the average F at equilibrium that all F's will be scaled to in order
      #to achieve equal F in every year. As depletion approaches the target value this 
      #F will approach F(OFL).
      F_report<-SPRfit$F_report
      FScale<-mean(F_report[(length(F_report)-9):length(F_report)])
      if(!is.null(ABC_Fraction)){
        F.ABC<-ABC_Fraction*FScale
      }else{
        F.ABC<-FScale
      }
      
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
        
        Depletion <- TimeFit3$SpawnBio/Virgin_bio
        Target.Rebuild <- mean(Depletion[(length(Depletion)-9):length(Depletion)])
      }else if(Forecast_target==3){
        Target.Depletion <- forecast[["Btarget"]]
        Depletion<-TimeFit3$SpawnBio/Virgin_bio

        Achieved.Depletion <- mean(Depletion[(length(Depletion)-9):length(Depletion)])
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
      }
    }else if(fitting_ABC==TRUE){
      DepletionScale<-1
      FScale<-F.ABC
    }else if(fitting_Rebuild==TRUE){
      DepletionScale<-1
      FScale<-F.ABC
      F_report<-SPRfit$F_report
      F_Rebuild_Scale<-F_report[SPRfit$Yr==Rebuild_yr]
      
      if(Forecast_target==1){
        Target.Rebuild <- forecast[["SPRtarget"]]
        Depletion<-SPRfit$SPR
        
        Achieved.Rebuild <- mean(Depletion[SPRfit$Yr==Rebuild_yr])
        Rebuild.Scale <- (1-Target.Rebuild)/(1-Achieved.Rebuild)
        Rebuild.Scale <- (-log(1-((1-exp(-F_Rebuild_Scale))*Rebuild.Scale))/F_Rebuild_Scale)
        Rebuild.Scale <- Rebuild.Scale*F_Rebuild_Scale
      }else if(Forecast_target==2){
        Depletion<-TimeFit3$SpawnBio/Virgin_bio
        
        Achieved.Rebuild <- mean(Depletion[SPRfit$Yr==Rebuild_yr])
        Rebuild.Scale <- (1-Target.Rebuild)/(1-Achieved.Rebuild)
        Rebuild.Scale <- (-log(1-((1-exp(-F_Rebuild_Scale))*Rebuild.Scale))/F_Rebuild_Scale)
        Rebuild.Scale <- Rebuild.Scale*F_Rebuild_Scale
      }else if(Forecast_target==3){
        Target.Rebuild <- forecast[["Btarget"]]
        Depletion<-TimeFit3$SpawnBio/Virgin_bio
        
        Achieved.Rebuild <- mean(Depletion[SPRfit$Yr==Rebuild_yr])
        Rebuild.Scale <- (1-Target.Rebuild)/(1-Achieved.Rebuild)
        Rebuild.Scale <- (-log(1-((1-exp(-F_Rebuild_Scale))*Rebuild.Scale))/F_Rebuild_Scale)
        Rebuild.Scale <- Rebuild.Scale*F_Rebuild_Scale
      }
      Rebuild.Scale <- min(Rebuild.Scale,FScale)
    }
    
    Fmult2 <- FScale/SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))]
    if(fitting_Rebuild==TRUE){
      Fmult2[rebuild_ref] <- Rebuild.Scale/SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))][rebuild_ref]
    }
    
    if(loop>1 & fitting_OFL==TRUE){
      F_adjust1 <- (Last_Mult1-1)/(Last_Mult1-DepletionScale)
      if(F_adjust1>0){
        if(DepletionScale<1){
          DepletionScale <- exp(log(DepletionScale)*F_adjust1)
        }
        if(DepletionScale>1){
          DepletionScale <- ((DepletionScale-1)*F_adjust1+1)
        }
      }
    }else if(loop>1 & fitting_ABC==TRUE){
      F_adjust2 <- (Last_Mult2-1)/(Last_Mult2-mean(Fmult2[-fixed_ref]))
      if(F_adjust2 > 0){
        if(length(Fmult2[Fmult2<1])>0){
          Fmult2[Fmult2<1] <- exp(log(Fmult2[Fmult2<1])*F_adjust2)
        }
        if(length(Fmult2[Fmult2>1])>0){
          Fmult2[Fmult2>1] <- ((Fmult2[Fmult2>1]-1)*F_adjust2+1)
        }
      }
    }else if(loop>1 & fitting_Rebuild==TRUE){
      F_adjust2a <- (Last_Mult2a-1)/(Last_Mult2a-mean(Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]]))
      F_adjust2b <- (Last_Mult2b-1)/(Last_Mult2b-mean(Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))]))
      if(!is.na(F_adjust2a)){
      if(F_adjust2a > 0){
        Fmult2a <- Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]]
        if(length(Fmult2a[Fmult2a<1])>0){
          Fmult2a[Fmult2a<1] <- exp(log(Fmult2a[Fmult2a<1])*F_adjust2a)
        }
        if(length(Fmult2a[Fmult2a>1])>0){
          Fmult2a[Fmult2a>1] <- ((Fmult2a[Fmult2a>1]-1)*F_adjust2a+1)
        }
        Fmult2[rebuild_ref[which(!is.element(rebuild_ref,fixed_ref))]] <- Fmult2a
      }}
      if(!is.na(F_adjust2b)){
      if(F_adjust2b > 0){
        Fmult2b <- Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))]
        if(length(Fmult2b[Fmult2b<1])>0){
          Fmult2b[Fmult2b<1] <- exp(log(Fmult2b[Fmult2b<1])*F_adjust2b)
        }
        if(length(Fmult2b[Fmult2b>1])>0){
          Fmult2b[Fmult2b>1] <- ((Fmult2b[Fmult2b>1]-1)*F_adjust2b+1)
        }
        Fmult2[-sort(unique(c(fixed_ref,rebuild_ref)))] <- Fmult2b
      }}
    }
    
    Fmult1 <- rep(DepletionScale,100*length(seasons)*length(F_cols))
    
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
    
    Fmult1[fixed_ref] <- 1
    Fmult2[fixed_ref] <- 1
    Fmult3[fixed_ref] <- 1
    Comb_Mult <- Fmult1*Fmult2*Fmult3
    
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
    
    if(max(abs(1-Fmult1))>=Depletion.Threshold | max(abs(1-Fmult2))>=Annual.F.Threshold | max(abs(1-Fmult3))>=Allocation.Threshold){keepFitting<-TRUE}else{keepFitting<-FALSE}
    if(fitting_OFL==TRUE & Forecast_target==2 & loop<10){keepFitting<-TRUE}
    
    forecast_F[,4] <- forecast_F[,4]*Comb_Mult
    forecast_F[fixed_ref,4] <- Fixed_catch_target[,4]
    forecast[["ForeCatch"]] <- forecast_F
    #Write the modified forecast data out to a file and rerun projections
    SS_writeforecast(mylist=forecast,overwrite = TRUE)
    shell(paste("cd /d ",getwd()," && ss -nohess",sep=""))
    #If all values have converged check if this is the OFL loop or the P* loop
    if(keepFitting==FALSE){
      if(fitting_OFL==TRUE){
        fitting_OFL <- FALSE
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
          #If in the P* loop then reset to keep fitting and change the target from P* to Rebuild
          loop <- 0
          method <- "Rebuild"
          keepFitting <- TRUE
          fitting_Rebuild <- TRUE
        }
        #If on the P* loop write out the P* results to a new folder (replace any old folder that exists)
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