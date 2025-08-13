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
#  Last update date: 8/1/25                                           #
#  Contact: nathan.vaughan@noaa.gov                                     #
#                                                                       #
##---------------------------------------------------------------------##
#
#   Disclaimer
#
#   "The United States Department of Commerce (DOC) GitHub project code 
#    is provided on an "as is" basis and the user assumes responsibility 
#    for its use. DOC has relinquished control of the information and no 
#    longer has responsibility to protect the integrity, confidentiality, 
#    or availability of the information. Any claims against the Department 
#    of Commerce stemming from the use of its GitHub project will be 
#    governed by all applicable Federal law. Any reference to specific 
#    commercial products, processes, or services by service mark, trademark, 
#    manufacturer, or otherwise, does not constitute or imply their endorsement, 
#    recommendation or favoring by the Department of Commerce. The Department 
#    of Commerce seal and logo, or the seal and logo of a DOC bureau, shall 
#    not be used in any manner to imply endorsement of any commercial 
#    product or activity by DOC or the United States Government."
#
#
##                                                                       #
##---------------------------------------------------------------------##
##

# This is now a generalized function to run projections that should work for most species. 
# The function is now designed to mostly read out settings from the forecast file the 
# same way you would generally set things up for a base SS projection.
# The only additional inputs needed are a proportion for F.ABC relative to F.OFL if 
# ABC projections are desired and a rebuild target year if rebuilding projections are desired.
#

run.projections<-function(Assessment_dir, #Here you set the location of a previously fit SS3.3 stock assessment to perform projections
                          ABC_Fraction = NULL, #Set the ABC target as a fraction of the OFL target if NULL will not fit ABC projections
                          Rebuild_yr = NULL, #Set the rebuild target year if NULL will not fit rebuild projections
                          Calc_F0 = FALSE, #Should an F=0 projection be performed
                          Const_Catch = NULL, #Constant catch target in mt. If constant catch is chosen no other simulations will run (i.e. benchmark, OFL, etc)  
                          F_max = FALSE, #If true and forecast method is Fmsy will replace Fmsy with Fmax search (maximum yield per recruit)
                          Depletion.Threshold = 0.0001, # These are all just thresholds for when 
                          Annual.F.Threshold = 0.0001, # targets are acceptably achieved these default to a .01% change
                          Allocation.Threshold = 0.0001, # increase them if run is too slow or reduce to improve fit if run is fast.
                          Step.Threshold = 0.0001, #
                          Benchmark_complete = FALSE, #Set to true if you already have a fit Benchmark run but want to add or adjust an ABC or Rebuild run
                          Benchmark_recruit_setting = 0, #Forecast recruitment setting to use for benchmarks in forecast file
                          OFL_recruit_setting = NULL, #Forecast recruitment setting to use for OFL in forecast file if NULL uses forecast file
                          ABC_recruit_setting = NULL, #Forecast recruitment setting to use for OFL in forecast file if NULL uses forecast file
                          Rebuild_recruit_setting = NULL, #Forecast recruitment setting to use for OFL in forecast file if NULL uses forecast file
                          F0_recruit_setting = NULL, #Forecast recruitment setting to use for OFL in forecast file if NULL uses forecast file
                          Const_catch_recruit_setting = NULL, #Forecast recruitment setting to use for OFL in forecast file if NULL uses forecast file
                          Make_plots = FALSE, #Should plots be created (this is useful for diagnostics but can cause annoying errors if plot window is small)
                          Calc_Hessian = FALSE, #Should the hessian inversion be completed for runs once converged. TODO: NOT YET IMPLEMENTED!!!
                          Do_Pstar = FALSE, #If TRUE then ABC_Fraction above will instead be the P* probability of overfishing limit for ABC calculation. TODO: NOT YET IMPLEMENTED!!!
                          Years_report = 20, #How many years of projection to include in stored OFL and ABC reporting. (All forecast years data will still be available in report file) 
                          Years_projection = 100, #How many years of projection to run (need enough to reach equilibrium) 100 is safe but may not be sufficient for some long lived species.
                          Constant_fixed_catch = NULL, #Input a data frame with column names ("Fleet","Catch or F","Basis")  
                          Annual_fixed_catch = NULL, #Input data frame for fixed catches in format of forecatch data frame with columns c("Year","Seas","Fleet","Catch or F","Basis")
                          Starting_Forecatch = NULL, #Input data frame for initial F values in format of forecatch data frame with columns c("Year","Seas","Fleet","Catch or F","Basis") otherwise will default to recent mean F
                          Fleet_group = NULL, #Data frame with columns c("Fleet","Group")specifying fleet grouping for allocations (defaults to forecast file settings if NULL)
                          Group_Allocations = NULL, #Data frame specifying scenarios for allocation of catch between groups (defaults to forecast file settings if NULL)
                          #If not null input a dataframe with a column for each group and row for each scenario.
                          Run_in_MSE = FALSE,
                          MSY_step = 0.1,
                          SS_exe = NULL,
                          Verbose = FALSE,
                          Messages = TRUE
                          ) 
{
  
  projection_results <- list()
  
  #SSMSE::report_message("Running the allocation forecasting function.") 
  
  #Removed these as inputs as they are not needed yet, could add back to input options later
  if(Messages == TRUE){
    message("Starting projections calculations")
  }
  oldwd<-getwd()
  return_warnings <- FALSE
  combinded_warnings <- NULL
  # library(r4ss)
  # if(Run_in_MSE==TRUE){
  #   library(SSMSE)
  # }
  #Set large max print to avoid issues with writing out a large forecast file.
  options(max.print = 1000000)
  setwd(normalizePath(paste0(Assessment_dir)))
  Assessment_dir <- normalizePath(getwd())
  #Read in all the model files and results 
  start <- r4ss::SS_readstarter(verbose = Verbose)
  dat <- r4ss::SS_readdat(file = start$datfile, version = 3.3, verbose = Verbose)
  ctl <- r4ss::SS_readctl(file = start$ctlfile, version = 3.3, use_datlist = TRUE, datlist = dat, verbose = Verbose)
  results <- r4ss::SS_output(dir = getwd(), covar = FALSE, verbose = Verbose, printstats= Verbose)
  forecast <- r4ss::SS_readforecast(verbose = Verbose) 
  
  if(!is.null(forecast$ForeCatch)){
    if(length(grep("year",colnames(forecast$ForeCatch)))==1){
    forecast$ForeCatch <- forecast$ForeCatch[forecast$ForeCatch$year<=(dat$endyr+forecast$Nforecastyrs) &
                                               forecast$ForeCatch$year>(dat$endyr),]
    }else if(length(grep("Year",colnames(forecast$ForeCatch)))==1){
      forecast$ForeCatch <- forecast$ForeCatch[forecast$ForeCatch$Year<=(dat$endyr+forecast$Nforecastyrs) &
                                                 forecast$ForeCatch$Year>(dat$endyr),]
    }else{
      stop("Error: for some reason the forecast ForeCatch dataframe doesn't 
            have any column named 'year' or 'Year' contact developer to update
            for an apparent change in SS/r4ss formating")
    }
  }
  
  os <- .Platform[["OS.type"]]
  if(Run_in_MSE == TRUE){
    bin <- SSMSE:::get_bin("ss")
  }else{
    if(!is.null(SS_exe)){
      
    }else if(file.exists("ss")){
      SS_exe <- "ss"
    }else if(file.exists("ss.exe")){
      SS_exe <- "ss"
    } else if(file.exists("ss3.exe")){
      SS_exe <- "ss3"
    } else if(file.exists("SS_opt.exe")){
      SS_exe <- "SS_opt"
    } else if(file.exists("ss3_opt.exe")){
      SS_exe <- "ss3_opt"
    } else if(file.exists("ss3_win.exe")){
      SS_exe <- "ss3_win"
    }else{
      stop("Error: Couldn't find an expected SS executable name")
    }
  }
  if(file.exists("ss.par")){
    par_name <- "ss.par"
  }else if(file.exists("ss3.par")){
    par_name <- "ss3.par"
  }else{
    stop("Error: No par file found with name ss.par or ss3.par")
  }
  admb_options <- "-nohess"
  parlist <- r4ss::SS_readpar_3.30(parfile = par_name, datsource = dat, ctlsource = ctl)
  
  if(!is.null(Const_Catch)){
    Const_Catch <- sort(Const_Catch)
  }
  
  #First set up a working director for running projections in (to avoid overwriting the base files with a failed model run)
  #then copy all of the assessment files to this working folder (ignore any output directories that have been previously created)
  if(dir.exists(file.path(getwd(),"Working_dir"))){
    unlink(file.path(getwd(),"Working_dir"), recursive = TRUE)
  }
  dir.create(file.path(getwd(),"Working_dir"))
  
  if(Benchmark_complete == FALSE){
  temp.files <- list.files(path=Assessment_dir)
  folders <- c(which(temp.files=="Benchmark_target"), which(temp.files=="OFL_target"), which(temp.files=="ABC_target"), which(temp.files=="Rebuild_target"), which(temp.files=="F0_target"), which(temp.files=="Working_dir"))
  if(length(folders)>0){
    temp.files <- temp.files[-folders]
  }
  file.copy(from = file.path(getwd(),temp.files), to = file.path(getwd(),"Working_dir",temp.files))
  }else{
    temp.files <- list.files(path=file.path(getwd(),"Benchmark_target"))
    folders <- which(grep("Allocation_run_",temp.files,fixed=TRUE))
    if(length(folders)>0){
      temp.files <- temp.files[-folders]
    }
    file.copy(from = file.path(getwd(),"Benchmark_target",temp.files), to = file.path(getwd(),"Working_dir",temp.files))
  }
  
  #Set the new working directory 
  setwd(file.path(getwd(),"Working_dir"))
  
  #Modify assessment files to produce expected timeseries length and outputs
  #need enough projection years to allow equilibrium to be achieved
  forecast[["Nforecastyrs"]] <- Years_projection
  
  #If fitting the Benchmark/OFL values then save the specified recruitment source for catches 
  #and set to 0 (S/R curve) in order to first calculate benchmarks
  input_recruit_setting <- forecast$fcast_rec_option
  if(is.null(Benchmark_recruit_setting)){
    Benchmark_recruit_setting <- input_recruit_setting
  }
  if(is.null(OFL_recruit_setting)){
    OFL_recruit_setting <- input_recruit_setting
  }
  if(is.null(ABC_recruit_setting)){
    ABC_recruit_setting <- input_recruit_setting
  }
  if(is.null(Rebuild_recruit_setting)){
    Rebuild_recruit_setting <- input_recruit_setting
  }
  if(is.null(F0_recruit_setting)){
    F0_recruit_setting <- input_recruit_setting
  }
  if(is.null(Const_catch_recruit_setting)){
    Const_catch_recruit_setting <- input_recruit_setting
  }
  
  if(Benchmark_complete==FALSE){
    forecast$fcast_rec_option <- Benchmark_recruit_setting
  }
  
  #Need to set recdevs and implementation error for all projection years 
  #so that reading from par file is possible
  expected_forecast_rec_length <- length((min(dat[["endyr"]],ctl[["MainRdevYrLast"]])+1):(dat[["endyr"]]+forecast[["Nforecastyrs"]]))
  if(!is.null(parlist$recdev_forecast)){
    if(length(parlist$recdev_forecast[,1])!=expected_forecast_rec_length){
      temp_recs<-parlist$recdev_forecast[,2]
      parlist$recdev_forecast <- matrix(NA, nrow = expected_forecast_rec_length, ncol = 2)
      parlist$recdev_forecast[,1] <- (min(dat[["endyr"]],ctl[["MainRdevYrLast"]])+1):(dat[["endyr"]]+forecast[["Nforecastyrs"]])
      parlist$recdev_forecast[,2] <- rep(0,expected_forecast_rec_length)
      parlist$recdev_forecast[1:length(temp_recs),2] <- temp_recs
      colnames(parlist$recdev_forecast) <- c("year","recdev")
    }
  }else{
    parlist$recdev_forecast <- matrix(NA, nrow = expected_forecast_rec_length, ncol = 2)
    parlist$recdev_forecast[,1] <- (min(dat[["endyr"]],ctl[["MainRdevYrLast"]])+1):(dat[["endyr"]]+forecast[["Nforecastyrs"]])
    parlist$recdev_forecast[,2] <- rep(0,expected_forecast_rec_length)
    colnames(parlist$recdev_forecast) <- c("year","recdev")
  }
  
  # if(!is.null(parlist$Fcast_impl_error)){
  parlist$Fcast_impl_error <- matrix(NA, nrow = forecast[["Nforecastyrs"]], ncol = 2)
  parlist$Fcast_impl_error[,1] <- (dat[["endyr"]]+1):(dat[["endyr"]]+forecast[["Nforecastyrs"]])
  parlist$Fcast_impl_error[,2] <- rep(0,forecast[["Nforecastyrs"]])
  colnames(parlist$Fcast_impl_error) <- c("year","impl_error")
  if(forecast[["stddev_of_log_catch_ratio"]]==0){
    forecast[["stddev_of_log_catch_ratio"]]<-0.001
  }
  
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
  Dead_cols <- grep("dead(B)", names(TimeFit), fixed = TRUE)
  CatchN_cols <- grep("retain(N)", names(TimeFit), fixed = TRUE)
  DeadN_cols <- grep("dead(N)", names(TimeFit), fixed = TRUE)
  F_cols <- grep("F", names(TimeFit), fixed = TRUE)
  if(length(F_cols)==0){
    F_cols <- grep("Hrate", names(TimeFit), fixed = TRUE)
  }
  TimeFit2 <- aggregate(TimeFit[,sort(c(2,4,7,Catch_cols,Dead_cols,CatchN_cols,DeadN_cols,F_cols))],by=list(TimeFit$Yr,TimeFit$Seas),FUN=sum,na.rm=TRUE)[,-c(3,4,5)]
  names(TimeFit2)[c(1,2)] <- c("Yr", "Seas")
  
  Catch_cols2 <- grep("retain(B)", names(TimeFit2), fixed = TRUE)
  Dead_cols2 <- grep("dead(B)", names(TimeFit2), fixed = TRUE)
  CatchN_cols2 <- grep("retain(N)", names(TimeFit2), fixed = TRUE)
  DeadN_cols2 <- grep("dead(N)", names(TimeFit2), fixed = TRUE)
  F_cols2 <- grep("F", names(TimeFit2), fixed = TRUE)
  if(length(F_cols2)==0){
    F_cols2 <- grep("Hrate", names(TimeFit2), fixed = TRUE)
  }
  
  TimeFit3 <- aggregate(TimeFit[,sort(c(2,7,8,Catch_cols,Dead_cols,CatchN_cols,DeadN_cols,F_cols))],by=list(TimeFit$Yr),FUN=sum,na.rm=TRUE)[,-2]
  names(TimeFit3)[c(1)] <- c("Yr")
  Virgin_bio <- TimeFit3$SpawnBio[1]
  
  Catch_cols3 <- grep("retain(B)", names(TimeFit3), fixed = TRUE)
  Dead_cols3 <- grep("dead(B)", names(TimeFit3), fixed = TRUE)
  CatchN_cols3 <- grep("retain(N)", names(TimeFit3), fixed = TRUE)
  DeadN_cols3 <- grep("dead(N)", names(TimeFit3), fixed = TRUE)
  F_cols3 <- grep("F", names(TimeFit3), fixed = TRUE)
  if(length(F_cols3)==0){
    F_cols3 <- grep("Hrate", names(TimeFit3), fixed = TRUE)
  }
  
  achieved.report <- TimeFit2[0,1:8]
  colnames(achieved.report)<-c("Year","Seas","Fleet","retain(B)","dead(B)","retain(N)","dead(N)","F")
  for(i in 1:length(F_cols2)){
    temp.data <- TimeFit2[,1:8]
    temp.data[,3] <- i 
    temp.data[,4] <- TimeFit2[,Catch_cols2[i]]
    temp.data[,5] <- TimeFit2[,Dead_cols2[i]]
    temp.data[,6] <- TimeFit2[,CatchN_cols2[i]]
    temp.data[,7] <- TimeFit2[,DeadN_cols2[i]]
    temp.data[,8] <- TimeFit2[,F_cols2[i]]
    colnames(temp.data)<-c("Year","Seas","Fleet","retain(B)","dead(B)","retain(N)","dead(N)","F")
    achieved.report <- rbind(achieved.report,temp.data)
  }
  achieved.report<-achieved.report[order(achieved.report$Year,achieved.report$Seas,achieved.report$Fleet),]
 
  #Set all future projections to fish at constant apical F that matches recent years
  #get the years of timeseries F's based on the forecast year range
  if(is.data.frame(forecast[["Fcast_years"]])){
    min_fcast_yr <- forecast[["Fcast_years"]][forecast[["Fcast_years"]][,'MG_type']==11,'st_year']
    if(min_fcast_yr==-999){
      min_fcast_yr <- dat[["styr"]]
    }else if(min_fcast_yr <= 0){
      min_fcast_yr <- dat[["endyr"]] + min_fcast_yr
    }
    
    max_fcast_yr <- forecast[["Fcast_years"]][forecast[["Fcast_years"]][,'MG_type']==11,'end_year']
    if(max_fcast_yr <= 0){
      max_fcast_yr <- dat[["endyr"]] + max_fcast_yr
    }
    
  }else{ 
    if(forecast[["Fcast_years"]][3]==-999){
      min_fcast_yr <- dat[["styr"]]
    }else if(forecast[["Fcast_years"]][3]>0){
      min_fcast_yr <- forecast[["Fcast_years"]][3]
    }else{
      min_fcast_yr <- dat[["endyr"]]+forecast[["Fcast_years"]][3]
    }
    
    if(forecast[["Fcast_years"]][4]>0){
      max_fcast_yr <- forecast[["Fcast_years"]][4]
    }else{
      max_fcast_yr <- dat[["endyr"]]+forecast[["Fcast_years"]][4]
    }
  }
  TargetYears <- TimeFit2[TimeFit2$Yr>=min_fcast_yr & TimeFit2$Yr<=max_fcast_yr,]
  TargetYears <- TargetYears[,c(2,F_cols2)]
  seasons <- unique(TargetYears[,1])
  F_by_Fleet_seas <- as.data.frame(matrix(apply(TargetYears[TargetYears[,1]==seasons[1],,drop=FALSE], 2, mean),nrow=1,ncol=(length(F_cols)+1)))
  if(length(seasons)>1){
    for(i in seasons[-1]){
      F_by_Fleet_seas <- rbind(F_by_Fleet_seas,apply(TargetYears[TargetYears[,1]==i,,drop=FALSE], 2, mean))
    }
  }
  Forecast_target <- forecast[["Forecast"]]
  if(!is.element(Forecast_target,c(1,2,3))){
    stop("forecast should be set to either 1, 2, or 3 so we know what the target is")
  }
  
  #Set up the forecast Forecatch dataframe to specify fixed fishing mortality rates
  #for each fleet for the entire 100 year projection period.
  #data frame will be build sequentially 
  #1) Set all F's to recent mean from the model
  Forecast_catch_setup<-matrix(1,nrow=forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols),ncol=7)
  Forecast_catch_setup[,1]<-sort(rep((dat[["endyr"]]+1):(dat[["endyr"]]+forecast[["Nforecastyrs"]]),length(seasons)*length(F_cols)))
  Forecast_catch_setup[,2]<-rep(sort(rep(seasons,length(F_cols))),forecast[["Nforecastyrs"]])
  Forecast_catch_setup[,3]<-rep(sort(which(dat$fleetinfo$type!=3)),forecast[["Nforecastyrs"]]*(length(seasons)))
  for(i in seasons){
    Forecast_catch_setup[Forecast_catch_setup[,2]==i,4]<-unlist(rep(F_by_Fleet_seas[F_by_Fleet_seas[,1]==i,-1],forecast[["Nforecastyrs"]]))
  }
  Forecast_catch_setup[,5] <- 99
  Forecast_catch_setup[,6] <- 0
  Forecast_catch_setup[,6] <- 0
  Forecast_catch_setup<-as.data.frame(Forecast_catch_setup)
  colnames(Forecast_catch_setup)<-c("Year","Seas","Fleet","Catch or F","Basis","Fixed","Rebuild")
  #2) Set any fleets with 0 F to a small number so they can be increased if needed
  Forecast_catch_setup[,"Catch or F"] <- ifelse(Forecast_catch_setup[,"Catch or F"]==0, 0.00000001,Forecast_catch_setup[,"Catch or F"])

  #3) Incorporate fixed values from the existing forecast file
  if(!is.null(forecast$ForeCatch)){
    for(i in seq_along(forecast$ForeCatch[,1])){
      match_row <- which(Forecast_catch_setup[,c("Year")]==forecast$ForeCatch[i,c("Year")] &
                         Forecast_catch_setup[,c("Seas")]==forecast$ForeCatch[i,c("Seas")] &
                         Forecast_catch_setup[,c("Fleet")]==forecast$ForeCatch[i,c("Fleet")])
      Forecast_catch_setup[match_row,"Catch or F"] <- forecast$ForeCatch[i,"Catch or F"]
      Forecast_catch_setup[match_row,"Fixed"] <- 1
      if(length(forecast$ForeCatch[i,])==5){
        Forecast_catch_setup[match_row,"Basis"] <- forecast$ForeCatch[i,"Basis"]
      }else if(length(forecast$ForeCatch[i,])==4){
        Forecast_catch_setup[match_row,"Basis"] <- forecast$InputBasis
      }else{
        stop("You have a Forecatch dataframe in your forecast file but it doesn't have 4 or 5 columns something is wrong")
      }
    }
    if(!is.null(Constant_fixed_catch)){
      return_warnings <- TRUE
      combinded_warnings <- paste(combinded_warnings,
                                  "You have input fixed constant catches but also have fixed values in the model forecast file. Check results to make sure all values are what you intended",
                                  sep="\n")
    }
    if(!is.null(Annual_fixed_catch)){
      return_warnings <- TRUE
      combinded_warnings <- paste(combinded_warnings,
                                  "You have input fixed annual catches but also have fixed values in the model forecast file. Check results to make sure all values are what you intended",
                                  sep="\n")
    }
  }
  #4) Fix constant F/Catch for fleets such as red tide, bycatch, or closed fleets  
  if(!is.null(Constant_fixed_catch)){
    for(i in seq_along(Constant_fixed_catch[,1])){
      Forecast_catch_setup[Forecast_catch_setup[,"Fleet"]==Constant_fixed_catch[i,"Fleet"],"Catch or F"] <- Constant_fixed_catch[i,"Catch or F"]
      Forecast_catch_setup[Forecast_catch_setup[,"Fleet"]==Constant_fixed_catch[i,"Fleet"],"Basis"] <- Constant_fixed_catch[i,"Basis"]
      Forecast_catch_setup[Forecast_catch_setup[,"Fleet"]==Constant_fixed_catch[i,"Fleet"],"Fixed"] <- 1
    }
  }
  #5) Update annual fixed F for fleets usually used for interim period catches
  if(!is.null(Annual_fixed_catch)){
    for(i in seq_along(Annual_fixed_catch[,1])){
      match_row <- which(Forecast_catch_setup[,c("Year")]==forecast$ForeCatch[i,c("Year")] &
                           Forecast_catch_setup[,c("Seas")]==forecast$ForeCatch[i,c("Seas")] &
                           Forecast_catch_setup[,c("Fleet")]==forecast$ForeCatch[i,c("Fleet")])
      Forecast_catch_setup[match_row,"Catch or F"] <- Annual_fixed_catch[i,"Catch or F"]
      Forecast_catch_setup[match_row,"Basis"] <- Annual_fixed_catch[i,"Basis"]
      Forecast_catch_setup[match_row,"Fixed"] <- 1
    }
  }
  #6) Create reference trackers for which years have fixed catch and which are adjusted
  fixed_ref <- which(Forecast_catch_setup[,"Fixed"]==1)
  adjusted_F_OFL<-which(Forecast_catch_setup[,"Fixed"]==0)
  Fixed_catch_target<-Forecast_catch_setup[fixed_ref,1:5]
  #7) Create forecast_F dataframe that will be used to update ForeCatch in forecast file
  forecast_F <- Forecast_catch_setup[,1:5]
  #6) Identify which years are subject to rebuilding plan F limits
  if(!is.null(Rebuild_yr)){
    #This reference identifies inputs subject to a rebuilding period F if any
    rebuild_ref <- which(Forecast_catch_setup[,1]<=Rebuild_yr)
    if(length(rebuild_ref)>0){
      Forecast_catch_setup[rebuild_ref,"Rebuild"]<-1
      adjusted_OFL_F_Rebuild<-adjusted_F_OFL[adjusted_F_OFL>max(rebuild_ref)]
      adjusted_Rebuild_F_Rebuild<-rebuild_ref
      if(!is.null(fixed_ref)){
        adjusted_Rebuild_F_Rebuild<-adjusted_Rebuild_F_Rebuild[-fixed_ref]
      }
    }
  }else{
    adjusted_OFL_F_Rebuild <- adjusted_F_OFL
    adjusted_Rebuild_F_Rebuild <- NULL
    rebuild_ref <- NULL
  }
  
  #Here all input values are assigned to an allocation group if needed and
  #relative landings targets are identified
  if(is.null(Fleet_group)){
    n_groups <- forecast[["N_allocation_groups"]]
  }else{
    n_groups <- length(unique(Fleet_group[Fleet_group[,"Group"]!=0,"Group"]))
  }
  groups <- rep(0,length(F_cols))
  fleets_by_group <- list()
  #Setup multiple allocation scenarios matrix to run allocation grid
  if(!is.null(Group_Allocations)){
    if(length(Group_Allocations[1,])!=n_groups){
      stop("Error: You have a different number of allocations specifed than the number of groups")
    }
    N_allocation_scenarios <- length(Group_Allocations[,1])
    Allocation_targets <- Group_Allocations
    for(i in seq_along(Allocation_targets[,1])){
      Allocation_targets[i,] <- Allocation_targets[i,]/sum(Allocation_targets[i,])
    }
    Allocation_targets <- rbind(unique(Fleet_group[Fleet_group[,"Group"]!=0,"Group"]),Group_Allocations)
  }else{
    N_allocation_scenarios <- 1
  }
  
  Allocation_tracker <- list()
  #Loop over all allocation scenarios
  for(i in 1:N_allocation_scenarios){
    projection_results[[paste0("Allocation_run_",i)]] <- list()
    Allocations <- forecast_F[,c(1,2,3,4,5,5)]
    Allocations[,c(4)] <- 0
    Allocations[,c(5,6)] <- 1
    names(Allocations) <- c("Year","Seas","Fleet","Group","Target","Achieved")
    
    if(n_groups>0){
      if(is.null(Fleet_group)){
        for(j in seq_along(forecast[["fleet_assignment_to_allocation_group"]][,"Fleet"])){
          groups[forecast[["fleet_assignment_to_allocation_group"]][j,"Fleet"]] <- forecast[["fleet_assignment_to_allocation_group"]][j,"Group"]
          Allocations[Allocations[,"Fleet"]==forecast[["fleet_assignment_to_allocation_group"]][j,"Fleet"],4] <- forecast[["fleet_assignment_to_allocation_group"]][j,"Group"]
        }
        alloc <- forecast[["allocation_among_groups"]][order(forecast[["allocation_among_groups"]][,"Year"]),]
        for(j in seq_along(alloc[,1])){
          for(k in seq_along(alloc[j,-1])){
            Allocations[Allocations[,1]>=alloc[,"Year"] & Allocations[,4]==k,5] <- alloc[j,(k+1)]/sum(alloc[j,-1])
          }
        }
      }else{
        for(j in seq_along(Fleet_group[,"Fleet"])){
          groups[Fleet_group[j,"Fleet"]] <- Fleet_group[j,"Group"]
          Allocations[Allocations[,"Fleet"]==Fleet_group[j,"Fleet"],4] <- Fleet_group[j,"Group"]
        }
        for(j in seq_along(Allocation_targets[1,])){
          Allocations[Allocations[,"Group"]==Allocation_targets[1,j],5] <- Allocation_targets[i+1,j]
        }
      }
      for(j in 1:n_groups){
        fleets_by_group[[j]]<-which(is.element(which(dat$fleetinfo$type!=3),which(groups==j)))
      }
    }
    Allocation_tracker[[i]] <- Allocations
  }
  Allocations <- Allocation_tracker[[1]]
  allocation_loop <- 1
 
  forecast[["Forecast"]] <- 4
  forecast[["InputBasis"]] <- -1
  forecast[["ForeCatch"]] <- forecast_F
  forecast[["FirstYear_for_caps_and_allocations"]] <- (dat[["endyr"]]+forecast[["Nforecastyrs"]]+1)
  forecast[["N_forecast_loops"]] <- 2
  
  last_forecast_F <- forecast[["ForeCatch"]]
  
  keepFitting <- TRUE
  loop <- 0
  subloop <- 0
  Curr_max_mult <- Last_max_mult <- Min_max_mult <- F_maxed <- 100000
  
  global_adjuster <- 1
  max_F_limit <- ctl$maxF
  F_adjust1 <- F_adjust1_2 <- 1
  F_adjust2 <- F_adjust2_2 <- F_SS_adjust <- F_adjust3 <- rep(1,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
  search_step <- MSY_step 
  Fmult1 <- Fmult2 <- Fmult3 <- Fmult4 <- rep(1.01,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
  Fmult1_raw <- Fmult2_raw <- Fmult3_raw <- Fmult4_raw <- rep(1.01,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
  Fmult2a <- Fmult2b <- 1
  First_run<-TRUE
  F_SS_adjust_year <- list(a=sort(rep(1:forecast[["Nforecastyrs"]],length(seasons)*length(F_cols))))
  
  if(!is.null(Const_Catch)){
    forecast$fcast_rec_option <- Const_catch_recruit_setting
    r4ss::SS_writepar_3.30(parlist = parlist,outfile=par_name,overwrite = TRUE, verbose = Verbose)
    r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
    r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
    
    dir <- normalizePath(getwd())
    if(Run_in_MSE==FALSE){
      bin <- file.path(dir,SS_exe)
    }
    if (os == "unix") {
      system(
        paste0(
          "cd ", dir, ";", paste0(bin, " "),
          admb_options
        ),
        ignore.stdout = TRUE
      )
    } else {
      system(paste0(paste0(bin, " "), admb_options),
             invisible = TRUE, ignore.stdout = TRUE,
             show.output.on.console = FALSE
      )
    }
    Sys.sleep(0.05)
    
    #Begin the search in the Benchmark phase
    fitting_Benchmark <- FALSE
    fitting_OFL <- FALSE
    fitting_ABC <- FALSE
    fitting_Rebuild <- FALSE
    fitting_F0 <- FALSE
    fitting_Fixed_Catch <- TRUE
    CC_loop <- 1
    Catch_Target <- rep(Const_Catch[CC_loop],forecast[["Nforecastyrs"]])
    Catch_trunc <- 0
    method <- "fixed_catch"
  }else if(Benchmark_complete == FALSE){
    #Save all the modified files and then perform a base run of SS so that output is specified correctly with 
    #a forecast[["Nforecastyrs"]] year projection series.
    
    r4ss::SS_writepar_3.30(parlist = parlist, outfile = par_name, overwrite = TRUE, verbose = Verbose)
    r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
    r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
    
    dir <- normalizePath(getwd())
    if(Run_in_MSE==FALSE){
      bin <- file.path(dir,SS_exe)
    }
    if (os == "unix") {
      system(
        paste0(
          "cd ", dir, ";", paste0(bin, " "),
          admb_options
        ),
        ignore.stdout = TRUE
      )
    } else {
      system(paste0(paste0(bin, " "), admb_options),
             invisible = TRUE, ignore.stdout = TRUE,
             show.output.on.console = FALSE
      )
    }
    Sys.sleep(0.05)
    
    
    #Begin the search in the Benchmark phase
    fitting_Benchmark <- TRUE
    fitting_OFL <- FALSE
    fitting_ABC <- FALSE
    fitting_Rebuild <- FALSE
    fitting_F0 <- FALSE
    fitting_Fixed_Catch <- FALSE
    
    MSY.Fit <- data.frame(catch=c(0),Ave.F=c(0),depletion=c(0),target.depletion=c(0))
    method <- "Benchmark"
    
  }else{
    #Save all the modified files and then set search to begin in OFL phase 
    
    start <- r4ss::SS_readstarter(verbose = Verbose)
    dat <- r4ss::SS_readdat(file = start$datfile, version = 3.3, verbose = Verbose)
    ctl <- r4ss::SS_readctl(file = start$ctlfile, version = 3.3, use_datlist = TRUE, datlist = dat, verbose = Verbose)
    results <- r4ss::SS_output(dir = getwd(), covar = FALSE, verbose = Verbose, printstats= Verbose)
    forecast <- r4ss::SS_readforecast(verbose = Verbose) 
    parlist <- r4ss::SS_readpar_3.30(parfile = par_name, datsource = dat, ctlsource = ctl, verbose = Verbose)
    
    forecast$fcast_rec_option <- OFL_recruit_setting
    r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
    
    fitting_Benchmark <- FALSE
    fitting_OFL <- TRUE
	  fitting_ABC <- FALSE
    fitting_Rebuild <- FALSE
    fitting_F0 <- FALSE
    fitting_Fixed_Catch <- FALSE
    method <- "OFL"
  }
  
  #Set up plot window for production of search diagnostic plots depending on target specifications
  #These plots were largely for diagnostic testing during code development but can allow you to see what
  #is going wrong if a future bug does occur.
  if(Make_plots==TRUE){
    if(Forecast_target==2 & fitting_Benchmark==TRUE){
      par(mfrow=c(5,2))
    }else{
      par(mfrow=c(4,2))
    }
  }
  
  if(Messages == TRUE){
    if(Benchmark_complete==FALSE){
      message("Projections targets have been set up begining optimization of Fs for Benchmarks")
    }else{
      message("Projections targets have been set up begining optimization of Fs for OFL")
    }
  }
  
  
  #Now start a loop of projecting and modifying fixed F's until the desired 
  #landings projections are achieved
  while(keepFitting){
    #Read in the SS results for landings and stock status to determine if desired
    #targets have been achieved
    resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
    TimeFit <- resultsFit[["timeseries"]]
    TimeFit <- TimeFit[TimeFit[,"Yr"]>dat[["endyr"]],]
    SPRfit <- resultsFit[["sprseries"]]
    SPRfit <- SPRfit[SPRfit[,"Yr"]>dat[["endyr"]],]
    loop <- loop + 1
    if(Make_plots==TRUE){
      try(
        {
          par(mar=c(4,3,3,2))
          plot(SPRfit[SPRfit[,"Yr"]>=dat[["endyr"]],"F_report"],xlab="year",ylab="F",main = paste0(method," loop = ",loop))
          plot(SPRfit[SPRfit[,"Yr"]>=dat[["endyr"]],"Deplete"],xlab="year",ylab="Depletion",main = paste0(method," loop = ",loop))
        }
      )
    }
    
    if(is.element(loop,c(1,2,3,4,5,10,20,30,40,50,100,200,300,400,500,1000))){
      if(Messages == TRUE){
        message(paste0("Running optimization loop ",loop))
        if(loop>50){
          message("Optimization seems to be taking a while I would suggest checking the working directory report and/or forcast files to ensure this is converging")
        }
        if(loop==1000){
          message("Something is probably wrong are you really sure this is converging?")
        }
      }
    }
    
    #Identify the column numbers for Catch, F, SSB, Recruits, etc
    Catch_cols <- grep("retain(B)", names(TimeFit), fixed = TRUE)
    Dead_cols <- grep("dead(B)", names(TimeFit), fixed = TRUE)
    CatchN_cols <- grep("retain(N)", names(TimeFit), fixed = TRUE)
    DeadN_cols <- grep("dead(N)", names(TimeFit), fixed = TRUE)
    F_cols <- grep("F", names(TimeFit), fixed = TRUE)
    if(length(F_cols)==0){
      F_cols <- grep("Hrate", names(TimeFit), fixed = TRUE)
    }
    
    TimeFit2 <- aggregate(TimeFit[,sort(c(2,4,7,Catch_cols,Dead_cols,CatchN_cols,DeadN_cols,F_cols))],by=list(TimeFit$Yr,TimeFit$Seas),FUN=sum,na.rm=TRUE)[,-c(3,4,5)]
    names(TimeFit2)[c(1,2)] <- c("Yr", "Seas")
    
    Catch_cols2 <- grep("retain(B)", names(TimeFit2), fixed = TRUE)
    Dead_cols2 <- grep("dead(B)", names(TimeFit2), fixed = TRUE)
    CatchN_cols2 <- grep("retain(N)", names(TimeFit2), fixed = TRUE)
    DeadN_cols2 <- grep("dead(N)", names(TimeFit2), fixed = TRUE)
    F_cols2 <- grep("F", names(TimeFit2), fixed = TRUE)
    if(length(F_cols2)==0){
      F_cols2 <- grep("Hrate", names(TimeFit2), fixed = TRUE)
    }
    
    TimeFit3 <- aggregate(TimeFit[,sort(c(2,7,8,Catch_cols,Dead_cols,CatchN_cols,DeadN_cols,F_cols))],by=list(TimeFit$Yr),FUN=sum,na.rm=TRUE)[,-2]
    names(TimeFit3)[c(1)] <- c("Yr")
    
    Catch_cols3 <- grep("retain(B)", names(TimeFit3), fixed = TRUE)
    Dead_cols3 <- grep("dead(B)", names(TimeFit3), fixed = TRUE)
    CatchN_cols3 <- grep("retain(N)", names(TimeFit3), fixed = TRUE)
    DeadN_cols3 <- grep("dead(N)", names(TimeFit3), fixed = TRUE)
    F_cols3 <- grep("F", names(TimeFit3), fixed = TRUE)
    if(length(F_cols3)==0){
      F_cols3 <- grep("Hrate", names(TimeFit3), fixed = TRUE)
    }
    
    achieved.report <- TimeFit2[0,1:8]
    colnames(achieved.report)<-c("Year","Seas","Fleet","retain(B)","dead(B)","retain(N)","dead(N)","F")
    for(i in 1:length(F_cols2)){
      temp.data <- TimeFit2[,1:8]
      temp.data[,3] <- i 
      temp.data[,4] <- TimeFit2[,Catch_cols2[i]]
      temp.data[,5] <- TimeFit2[,Dead_cols2[i]]
      temp.data[,6] <- TimeFit2[,CatchN_cols2[i]]
      temp.data[,7] <- TimeFit2[,DeadN_cols2[i]]
      temp.data[,8] <- TimeFit2[,F_cols2[i]]
      colnames(temp.data)<-c("Year","Seas","Fleet","retain(B)","dead(B)","retain(N)","dead(N)","F")
      achieved.report <- rbind(achieved.report,temp.data)
    }
    achieved.report<-achieved.report[order(achieved.report$Year,achieved.report$Seas,achieved.report$Fleet),]
    
    if(fitting_Rebuild==TRUE){
      terminal_year <- min(length(SPRfit$Yr[SPRfit$Yr<=Rebuild_yr])+5, Years_projection) #min to ensure report years beyond projection length not requested
    }else{
      terminal_year <- min(Years_report, Years_projection) #min to ensure report years beyond projection length not requested
    }
    
    Achieved.Catch <- apply(TimeFit3[,Catch_cols3,drop=FALSE],1,sum)[1:terminal_year]
    Achieved.Catch.All <- apply(TimeFit3[,Catch_cols3,drop=FALSE],1,sum)
    Achieved.SSBratio <- TimeFit3$SpawnBio[1:terminal_year]/Virgin_bio
    Achieved.SPR <- TimeFit3$SpawnBio[1:terminal_year]/TimeFit3$Recruit_0[1:terminal_year]
    Achieved.SSB <- TimeFit3$SpawnBio[1:terminal_year]
    Achieved.Rec <- TimeFit3$Recruit_0[1:terminal_year]
    Achieved.F <- SPRfit$F_report[1:terminal_year]
    if(is.null(Achieved.F)){
      Achieved.F <- SPRfit$F_std[1:terminal_year]
    }
    #If reading in results of a previous OFL run then set the F_OFL and F.ABC values before begining ABC/Rebuild loops
   
    if(Benchmark_complete==TRUE & First_run==TRUE){
      F_report<-SPRfit$F_report
      if(is.null(F_report)){
        F_report<-SPRfit$F_std
      }
      FScale<-median(F_report[(length(F_report)-0.5*Years_projection):length(F_report)])
      F_OFL<-FScale
      if(!is.null(ABC_Fraction)){
        F.ABC<-ABC_Fraction*FScale
      }else{
        F.ABC<-FScale
      }
      First_run<-FALSE
      Depletion<-TimeFit3$SpawnBio/Virgin_bio
      Target.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
      Target.Rebuild <- median(Depletion[(length(Depletion)-29):length(Depletion)])
      
      Achieved.Catch.equil <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/10
      
      if(n_groups>1){
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Benchmark"]]<-list()
        for(i in 1:n_groups){
          Achieved.Catch.group.equil <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3[fleets_by_group[[i]]]])/10
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Equil_Benchmark"]][[i]]<-Achieved.Catch.group.equil
          
          Achieved.Catch.group <- apply(TimeFit3[1:terminal_year,Catch_cols3[fleets_by_group[[i]]],drop=FALSE],1,sum)
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Benchmark"]][[i]]<-Achieved.Catch.group
        }
      }
      
      Achieved.SSBratio.equil <- median(TimeFit3$SpawnBio[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)])/Virgin_bio
      Achieved.SPR.equil <- median(TimeFit3$SpawnBio[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)]/TimeFit3$Recruit_0[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)]) #median(SPRfit$SPR[(length(SPRfit$SPR)-9):length(SPRfit$SPR)])
      Achieved.SSB.equil <- median(TimeFit3$SpawnBio[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)])
      Achieved.Rec.equil <- median(TimeFit3$Recruit_0[(length(TimeFit3$Recruit_0)-9):length(TimeFit3$Recruit_0)])
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_equil_Benchmark"]] <- Achieved.Catch.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_equil_Benchmark"]] <- F_OFL
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_equil_Benchmark"]] <- Achieved.SSBratio.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_equil_Benchmark"]] <- Achieved.SSB.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_equil_Benchmark"]] <- Achieved.SPR.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_equil_Benchmark"]] <- Achieved.Rec.equil
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Forecatch_Benchmark"]] <- achieved.report
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_Benchmark"]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_Benchmark"]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_Benchmark"]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_Benchmark"]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_Benchmark"]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_Benchmark"]] <- Achieved.Rec
    }
    
    if(is.na(max(abs(achieved.report[,'F']-forecast_F[,"Catch or F"])[adjusted_F_OFL]))){
      if(fitting_Fixed_Catch==TRUE){
        if(Catch_trunc >= (forecast$Nforecastyrs - 20)){
          Catch_trunc <- Catch_trunc + 1
        }else{
          Catch_trunc <- Catch_trunc + 5
        }
        Catch_Target[(forecast$Nforecastyrs-c((Catch_trunc-1):0))] <- 0
        forecast_F[(length(forecast_F[,1])-c((Catch_trunc*length(seasons)*length(F_cols)-1):0)),4] <- 0
      }
      forecast_F[,4] <- forecast_F[,4]*0.5
      loop <- loop - 1
    }else{
      
    if(max(abs(achieved.report[,'F']-forecast_F[,"Catch or F"])[adjusted_F_OFL])>0.1){
       if(fitting_Fixed_Catch==TRUE){
         if(Catch_trunc >= (forecast$Nforecastyrs - 20)){
           Catch_trunc <- Catch_trunc + 1
         }else{
           Catch_trunc <- Catch_trunc + 5
         }
         Catch_Target[(forecast$Nforecastyrs-c((Catch_trunc-1):0))] <- 0
         forecast_F[(length(forecast_F[,1])-c((Catch_trunc*length(seasons)*length(F_cols)-1):0)),4] <- 0
       }
       expected_annual <- 1-exp(-aggregate(forecast_F[,"Catch or F"],by=F_SS_adjust_year,FUN=sum)$x)
       achieved_annual <- 1-exp(-aggregate(achieved.report[,'F'],by=F_SS_adjust_year,FUN=sum)$x)
       F_SS_adjust <- rep(achieved_annual/expected_annual,each=length(seasons)*length(F_cols))
       F_maxed <- max(achieved.report[,'F'])
    }else{
      F_SS_adjust <- rep(1,length(forecast_F[,1]))
    }
    {
    if(fitting_Benchmark==TRUE){
      
      #Calculate the average F at equilibrium that all F's will be scaled to in order
      #to achieve equal F in every year. As depletion approaches the target value this 
      #F will approach F(OFL).
       
      
      F_report<-SPRfit$F_report
      if(is.null(F_report)){
        F_report<-SPRfit$F_std
      }
      FScale<-median(F_report[(length(F_report)-0.5*Years_projection):length(F_report)])
      F_OFL<-FScale
      if(!is.null(ABC_Fraction)){
        F.ABC<-ABC_Fraction*FScale
      }else{
        F.ABC<-FScale
      }
      
      Achieved.Catch.equil <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/10
      Achieved.SSBratio.equil <- median(TimeFit3$SpawnBio[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)])/Virgin_bio
      Achieved.SPR.equil <- median(TimeFit3$SpawnBio[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)]/TimeFit3$Recruit_0[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)]) #median(SPRfit$SPR[(length(SPRfit$SPR)-9):length(SPRfit$SPR)])
      Achieved.SSB.equil <- median(TimeFit3$SpawnBio[(length(TimeFit3$SpawnBio)-9):length(TimeFit3$SpawnBio)])
      Achieved.Rec.equil <- median(TimeFit3$Recruit_0[(length(TimeFit3$Recruit_0)-9):length(TimeFit3$Recruit_0)])
      
      
      
      if(n_groups>1){
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Benchmark"]]<-list()
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Equil_Benchmark"]]<-list()
        for(i in 1:n_groups){
          Achieved.Catch.group.equil <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3[fleets_by_group[[i]]]])/10
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Equil_Benchmark"]][[i]]<-Achieved.Catch.group.equil
          
          Achieved.Catch.group <- apply(TimeFit3[1:terminal_year,Catch_cols3[fleets_by_group[[i]]],drop=FALSE],1,sum)
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Benchmark"]][[i]]<-Achieved.Catch.group
        }
      }
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_equil_Benchmark"]] <- Achieved.Catch.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_equil_Benchmark"]] <- F_OFL
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_equil_Benchmark"]] <- Achieved.SSBratio.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_equil_Benchmark"]] <- Achieved.SSB.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_equil_Benchmark"]] <- Achieved.SPR.equil
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_equil_Benchmark"]] <- Achieved.Rec.equil
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Forecatch_Benchmark"]] <- achieved.report
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_Benchmark"]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_Benchmark"]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_Benchmark"]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_Benchmark"]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_Benchmark"]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_Benchmark"]] <- Achieved.Rec
      #Calculate depletion target adjustment scale depending on the specified target (SPR ratio, SSB ratio, or true MSY)   
      if(Forecast_target==1){
        search_step<-0.00001
        Target.Depletion <- forecast[["SPRtarget"]]
        Depletion<-SPRfit$SPR
        
        Achieved.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
        
        Achieved.Depletion <- min(Achieved.Depletion,max(Target.Depletion,0.9))
        
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        
        if(FScale==0){
          if(DepletionScale<=1.0001){
            DepletionScale<-1
          }else{
            FScale<-0.0001
            F_OFL<-FScale
            if(!is.null(ABC_Fraction)){
              F.ABC<-ABC_Fraction*FScale
            }else{
              F.ABC<-FScale
            }
          }
        }else{
          DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
        }
        
        Depletion_R<-TimeFit3$SpawnBio/Virgin_bio
        Target.Rebuild <- median(Depletion_R[(length(Depletion_R)-9):length(Depletion_R)])
        
      }else if(Forecast_target==2){
        Depletion <- TimeFit3$SpawnBio/Virgin_bio
        Achieved.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
        Achieved.Depletion <- min(Achieved.Depletion,.9)
        if(First_run == TRUE){
          Target.Depletion <- forecast[["Btarget"]]
          First_run <- FALSE
        }
        Target.Rebuild <- Target.Depletion
        Achieved.SSB <- Achieved.Depletion
        if(max(abs(1-Fmult3_raw))>Allocation.Threshold | 
		       max(abs(1-Fmult2_raw))>Annual.F.Threshold | 
		       max(abs(1-Fmult1_raw))>Depletion.Threshold){
          loop<-loop-1
          subloop<-subloop+1
          if(F_max==TRUE){
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/
			                  sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),"Recruit_0"])
          }else{
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/10
          }
          MSY.Fit[1,] <- c(Achieved.Catch,FScale,Achieved.Depletion,Target.Depletion)
        }else{
          subloop<-0
          if(F_max==TRUE){
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/
			                  sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),"Recruit_0"])
          }else{
            Achieved.Catch <- sum(TimeFit3[(length(TimeFit3[,1])-9):length(TimeFit3[,1]),Catch_cols3])/10
          }
          MSY.Fit <- rbind(MSY.Fit[1,],MSY.Fit)
          MSY.Fit[1,] <- c(Achieved.Catch,FScale,Achieved.Depletion,Target.Depletion)
          if(loop>1){
            if(Achieved.Catch<Last_Achieved_Catch){
              search_step <- -0.3*search_step
            }
          
            Target.Depletion <- Target.Depletion+search_step
            
            min_diff <- which(abs(MSY.Fit[,4]-Target.Depletion)<0.001)
            if(length(min_diff)>0){
              Old.Catch <- MSY.Fit[min_diff[1],1]
              if(Old.Catch<Achieved.Catch){
                search_step <- -0.3*search_step
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
        
        if(FScale==0){
          if(DepletionScale<=1.0001){
            DepletionScale<-1
          }else{
            FScale<-0.0001
            F_OFL<-FScale
            if(!is.null(ABC_Fraction)){
              F.ABC<-ABC_Fraction*FScale
            }else{
              F.ABC<-FScale
            }
          }
        }
        
        if(Make_plots==TRUE){
          try(
            {
              par(mar=c(4,3,3,2))
              if(F_max==TRUE){
                plot(x=TimeFit3[,"Yr"],y=apply(TimeFit3[,Catch_cols3,drop=FALSE],1,sum)/TimeFit3[,"Recruit_0"],
      		       xlab="year",ylab="Total Yield Per Recruit",main = paste0(method," loop = ",loop,".",subloop))
                plot(x=MSY.Fit[,"depletion"],y=MSY.Fit[,"catch"],xlim=c(0.9*min(MSY.Fit[,3:4]),1.1*max(MSY.Fit[,3:4])),
      		       xlab="Esimated depletion",ylab="Total Yield Per Recruit",main = paste0(method," loop = ",loop,".",subloop))
                lines(x=c(MSY.Fit[1,c(4,4)]),y=c(0,2*max(MSY.Fit[,"catch"])),col="dark red")
                points(x=MSY.Fit[1,3],y=MSY.Fit[1,1],pch=16,col="dark blue")
              }else{
                plot(x=TimeFit3[,"Yr"],y=apply(TimeFit3[,Catch_cols3,drop=FALSE],1,sum),
      		       xlab="year",ylab="Total Yield",main = paste0(method," loop = ",loop,".",subloop))
                plot(x=MSY.Fit[,"depletion"],y=MSY.Fit[,"catch"],xlim=c(0.9*min(MSY.Fit[,3:4]),1.1*max(MSY.Fit[,3:4])),
      		       xlab="Esimated depletion",ylab="Total Yield",main = paste0(method," loop = ",loop,".",subloop))
                lines(x=c(MSY.Fit[1,c(4,4)]),y=c(0,2*max(MSY.Fit[,"catch"])),col="dark red")
                points(x=MSY.Fit[1,3],y=MSY.Fit[1,1],pch=16,col="dark blue")
              }
            }
          )
        }
      }else if(Forecast_target==3){
        search_step<-0.00001
        Target.Depletion <- forecast[["Btarget"]]
        Target.Rebuild <- forecast[["Btarget"]]
        Depletion<-TimeFit3$SpawnBio/Virgin_bio
        Achieved.Depletion <- median(Depletion[(length(Depletion)-29):length(Depletion)])
        Achieved.Depletion <- min(Achieved.Depletion,max(0.9,Target.Depletion))
        DepletionScale <- (1-Target.Depletion)/(1-Achieved.Depletion)
        
        if(FScale==0){
          if(DepletionScale<=1.0001){
            DepletionScale<-1
          }else{
            FScale<-0.0001
            F_OFL<-FScale
            if(!is.null(ABC_Fraction)){
              F.ABC<-ABC_Fraction*FScale
            }else{
              F.ABC<-FScale
            }
          }
        }else{
          DepletionScale <- (-log(1-((1-exp(-FScale))*DepletionScale))/FScale)
        }
        
      }
    }else if(fitting_OFL==TRUE){
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the Benchmark MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that Benchmark search is complete
      FScale<-F_OFL #Set the F target to F.OFL for rescaling annual F values
      
      if(n_groups>1){
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_OFL"]]<-list()
        for(i in 1:n_groups){
          Achieved.Catch.group <- apply(TimeFit3[1:terminal_year,Catch_cols3[fleets_by_group[[i]]],drop=FALSE],1,sum)
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_OFL"]][[i]]<-Achieved.Catch.group
        }
      }
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_OFL"]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_OFL"]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_OFL"]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_OFL"]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_OFL"]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_OFL"]] <- Achieved.Rec
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Forecatch_OFL"]] <- achieved.report
      
    }else if(fitting_ABC==TRUE){
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the Benchmark MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that Benchmark search is complete
      FScale<-F.ABC #Set the F target to F.ABC for rescaling annual F values
      
      if(n_groups>1){
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_ABC"]]<-list()
        for(i in 1:n_groups){
          Achieved.Catch.group <- apply(TimeFit3[1:terminal_year,Catch_cols3[fleets_by_group[[i]]],drop=FALSE],1,sum)
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_ABC"]][[i]]<-Achieved.Catch.group
        }
      }
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_ABC"]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_ABC"]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_ABC"]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_ABC"]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_ABC"]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_ABC"]] <- Achieved.Rec
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Forecatch_ABC"]] <- achieved.report
    }else if(fitting_F0==TRUE){
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the Benchmark MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that Benchmark search is complete
      FScale<-0 #Set the F target to 0 for rescaling annual F values
      
      if(n_groups>1){
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_F0"]]<-list()
        for(i in 1:n_groups){
          Achieved.Catch.group <- apply(TimeFit3[1:terminal_year,Catch_cols3[fleets_by_group[[i]]],drop=FALSE],1,sum)
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_F0"]][[i]]<-Achieved.Catch.group
        }
      }
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_F0"]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_F0"]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_F0"]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_F0"]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_F0"]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_F0"]] <- Achieved.Rec
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Forecatch_F0"]] <- achieved.report
    }else if(fitting_Rebuild==TRUE){
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the Benchmark MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that Benchmark search is complete
      FScale<-F_OFL #Set the F target to F_OFL for rescaling annual F values in years after the rebuild period.
      F_report<-SPRfit$F_report
      if(is.null(F_report)){
        F_report<-SPRfit$F_std
      }
      F_Rebuild_Scale<-F_report[SPRfit$Yr==Rebuild_yr]
      Depletion<-TimeFit3$SpawnBio/Virgin_bio
      Achieved.Rebuild <- mean(Depletion[SPRfit$Yr==Rebuild_yr])
      
      Rebuild.Scale <- (1-Target.Rebuild)/(1-Achieved.Rebuild)
      Rebuild.Ratio <- Rebuild.Scale
      Rebuild.Scale <- min(-log(1-((1-exp(-F_Rebuild_Scale))*Rebuild.Scale)),FScale)
      if(Rebuild.Scale < 0.00001 & Rebuild.Ratio < 1){
        Rebuild.Scale <- 0
      }
      if(n_groups>1){
        projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Rebuild"]]<-list()
        for(i in 1:n_groups){
          Achieved.Catch.group <- apply(TimeFit3[1:terminal_year,Catch_cols3[fleets_by_group[[i]]],drop=FALSE],1,sum)
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Group_Catch_Rebuild"]][[i]]<-Achieved.Catch.group
        }
      }
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Catch_Rebuild"]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_Rebuild"]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Depletion_Rebuild"]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_Rebuild"]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["SPR_Rebuild"]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recruitment_Rebuild"]] <- Achieved.Rec
      projection_results[[paste0("Allocation_run_",allocation_loop)]][["Forecatch_Rebuild"]] <- achieved.report
    }else if (fitting_Fixed_Catch==TRUE){
      
      search_step<-0.00001 #Set search step to small value so it doesn't trigger continued loops this value is only needed during the Benchmark MSY search
      DepletionScale<-1 #Set depletion scale to 1 so it doesn't trigger continued loops now that Benchmark search is complete
      FScale<-0 #Set the F target to F_OFL for rescaling annual F values in years after the rebuild period.
      
      Fmult4 <- rep(Catch_Target/Achieved.Catch.All,each=(length(seasons)*length(F_cols)))
      
      Fmult4 <- ifelse(forecast_F[,4]>=max_F_limit,ifelse(Fmult4>1,1,Fmult4),Fmult4)
      
      Fmult4[is.na(Fmult4)] <- 1
      
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("Catch_FixedCatch_",CC_loop)]] <- Achieved.Catch
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("F_FixedCatch_",CC_loop)]] <- Achieved.F
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("Depletion_FixedCatch_",CC_loop)]] <- Achieved.SSBratio
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("SSB_FixedCatch_",CC_loop)]] <- Achieved.SSB
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("SPR_FixedCatch_",CC_loop)]] <- Achieved.SPR
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("Recruitment_FixedCatch_",CC_loop)]] <- Achieved.Rec
      projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("Forecatch_FixedCatch_",CC_loop)]] <- achieved.report
    }
      if(is.infinite(DepletionScale)){
        DepletionScale <- 1
      }
      if(is.na(DepletionScale)){
        DepletionScale <- 1
      }
      if(is.nan(DepletionScale)){
        DepletionScale <- 1
      }
      if(DepletionScale <= 0){
        DepletionScale <- 1
      }
      if(DepletionScale <= 0.5){
        DepletionScale <- 0.5
      }
      if(DepletionScale >= 2){
        DepletionScale <- 2
      }  
    Fmult1_raw <- rep(DepletionScale,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
    #Fmult2 calculations define the multiplier for adjusting annual F values
    #Zero catch years are identified first to prevent divide by zero errors in the scaling and
    #to tell the search algorithm that the target has been achieved

    if(is.null(SPRfit$F_report)){
      zero_catch <- which(SPRfit$F_std[sort(rep(seq_along(SPRfit$F_std),length(seasons)*length(F_cols)))]==0)
    }else{
      zero_catch <- which(SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))]==0)
    }
    if(fitting_Fixed_Catch==FALSE){
      if(length(zero_catch)>0){
        if(FScale==0){
          Fmult2[zero_catch] <- 1
          Fmult2[-zero_catch] <- 0
        }else{
          Fmult2[zero_catch] <- 2
          if(is.null(SPRfit$F_report)){
            temp_F <- SPRfit$F_std[sort(rep(seq_along(SPRfit$F_std),length(seasons)*length(F_cols)))][-zero_catch]
          }else{
            temp_F <- SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))][-zero_catch]
          }
          temp_F[temp_F>max_F_limit] <- max_F_limit
          temp_F[temp_F<(min(0.001,FScale))] <- min(0.001,FScale)
          Fmult2[-zero_catch] <- FScale/temp_F
        }
      }else if(FScale==0){
        Fmult2[] <- 0
      }else{
        if(is.null(SPRfit$F_report)){
          temp_F <- SPRfit$F_std[sort(rep(seq_along(SPRfit$F_std),length(seasons)*length(F_cols)))]
        }else{
          temp_F <- SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))]
        }
        temp_F[temp_F>max_F_limit] <- max_F_limit
        temp_F[temp_F<(min(0.001,FScale))] <- min(0.001,FScale)
        Fmult2 <- FScale/temp_F
      }
    }else{
      Fmult2 <- rep(1,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
    }
    
    #If in a rebuild search phase the rebuild years are now adjusted independently of the later F_OFL years
    if(fitting_Rebuild==TRUE){
      if(is.null(SPRfit$F_report)){
        temp_F <- SPRfit$F_std[sort(rep(seq_along(SPRfit$F_std),length(seasons)*length(F_cols)))][adjusted_Rebuild_F_Rebuild]
      }else{
        temp_F <- SPRfit$F_report[sort(rep(seq_along(SPRfit$F_report),length(seasons)*length(F_cols)))][adjusted_Rebuild_F_Rebuild]
      }
      zero_rebuild <- which(temp_F==0)
      if(length(zero_rebuild)>0){
        if(Rebuild.Scale==0){
          Fmult2[adjusted_Rebuild_F_Rebuild[zero_rebuild]] <- 1
          Fmult2[adjusted_Rebuild_F_Rebuild[-zero_rebuild]] <- 0
        }else{
          Fmult2[adjusted_Rebuild_F_Rebuild[zero_rebuild]] <- 2
          temp_F[temp_F>max_F_limit] <- max_F_limit
          temp_F[temp_F<(min(0.001,Rebuild.Scale))] <- min(0.001,Rebuild.Scale)
          Fmult2[adjusted_Rebuild_F_Rebuild[-zero_rebuild]] <- Rebuild.Scale/temp_F[-zero_rebuild]
        }
      }else if(Rebuild.Scale==0){
        Fmult2[adjusted_Rebuild_F_Rebuild] <- 0
      }else{
        temp_F[temp_F>max_F_limit] <- max_F_limit
        temp_F[temp_F<(min(0.001,Rebuild.Scale))] <- min(0.001,Rebuild.Scale)
        Fmult2[adjusted_Rebuild_F_Rebuild] <- Rebuild.Scale/temp_F
      }
    }
    Fmult2[which(is.infinite(Fmult2))] <- 1
    Fmult2[which(is.nan(Fmult2))] <- 1
    Fmult2[which(is.na(Fmult2))] <- 1
    Fmult2[which(Fmult2<=0)] <- 1
    Fmult2[which(Fmult2>=2)] <- 2
    Fmult2[which(Fmult2<=0.5)] <- 0.5
    Fmult2_raw <- Fmult2
    #Here a range of adjustments are made to the F step sizes based on the expected vs achieved change in F from
    #the previous step. i.e. if the last change only had half the impact expected on F then the next step will
    #be modified to make the next change in F twice as large as the raw change in F estimated.
    
    if(fitting_Fixed_Catch==FALSE){
      if((loop>1 |  subloop>2) & DepletionScale>0 & DepletionScale!=1 & fitting_Benchmark==TRUE & (Forecast_target!=2 | subloop>2)){
        F_adjust1 <- (4*F_adjust1 + global_adjuster)/5
        F_adjust1 <- F_adjust1*(Last_Mult1-1)/(Last_Mult1-DepletionScale)
        if(is.infinite(F_adjust1)){
          F_adjust1 <- global_adjuster
        }
        if(is.na(F_adjust1)){
          F_adjust1 <- global_adjuster
        }
        if(is.nan(F_adjust1)){
          F_adjust1 <- global_adjuster
        }
        if(F_adjust1 <= 0){
          F_adjust1 <- 0.5
        }
        if(F_adjust1 <= 0.5){
          F_adjust1 <- 0.5
        }
        if(F_adjust1 >= 2){
          F_adjust1 <- 2
        }
        if((DepletionScale-1)*(Last_Mult1-1) <= 0){
          F_adjust1_2 <- 0.8*F_adjust1_2
        }else{
          F_adjust1_2 <- (4*F_adjust1_2 + global_adjuster)/5
        }
        F_adjust1 <- F_adjust1 * F_adjust1_2
        if(F_adjust1>0){
          DepletionScale <- ((DepletionScale-1)*F_adjust1+1)
        }
        if(is.infinite(DepletionScale)){
          DepletionScale <- 1
        }
        if(is.na(DepletionScale)){
          DepletionScale <- 1
        }
        if(is.nan(DepletionScale)){
          DepletionScale <- 1
        }
        if(DepletionScale <= 0){
          DepletionScale <- 1
        }
        if(DepletionScale <= 0.5){
          DepletionScale <- 0.5
        }
        if(DepletionScale >= 2){
          DepletionScale <- 2
        }  
      }
      if(loop>1 & (fitting_Benchmark==TRUE | fitting_OFL==TRUE | fitting_ABC==TRUE)){
        F_adjust2 <- (4*F_adjust2 + global_adjuster)/5
        F_adjust2 <- F_adjust2*(Last_Mult2-1)/(Last_Mult2-F_SS_adjust*Fmult2)
        
        F_adjust2[which(is.infinite(F_adjust2))] <- global_adjuster
        F_adjust2[which(is.nan(F_adjust2))] <- global_adjuster
        F_adjust2[which(is.na(F_adjust2))] <- global_adjuster
        F_adjust2[which(F_adjust2<=0)] <- min(0.5,global_adjuster)
        F_adjust2[which(F_adjust2>=2)] <- 2
        #F_adjust2[which(F_adjust2<=0.5)] <- 0.5
        F_adjust2_2[which((Fmult2-1)*(Last_Mult2-1)<=0)] <- 0.8*F_adjust2_2[which((Fmult2-1)*(Last_Mult2-1)<=0)]
        #F_adjust2_2[which((Fmult2-1)*(Last_Mult2-1)>0)] <- 
        #  (4*F_adjust2_2[which((Fmult2-1)*(Last_Mult2-1)>0)]+global_adjuster)/5
        F_adjust2_2[which(abs(F_SS_adjust-1)>=0.1)] <- 0.8*F_adjust2_2[which(abs(F_SS_adjust-1)>=0.1)]
        
        Fmult2 <- ((Fmult2-1)*F_adjust2*F_adjust2_2+1)
        
        Fmult2[which(is.infinite(Fmult2))] <- 1
        Fmult2[which(is.nan(Fmult2))] <- 1
        Fmult2[which(is.na(Fmult2))] <- 1
        Fmult2[which(Fmult2<=0)] <- 1
        Fmult2[which(Fmult2>=2)] <- 2
        Fmult2[which(Fmult2<=0.5)] <- 0.5
      }
      if(loop>1 & min(Fmult2)>0 & fitting_Rebuild==TRUE & median(Fmult2[adjusted_Rebuild_F_Rebuild])!=1 & median(Fmult2[adjusted_OFL_F_Rebuild])!=1){
        if(exists("F_adjust2a")){
          F_adjust2a <- (F_adjust2a + global_adjuster)/2
        }else{
          F_adjust2a <- global_adjuster
        }
        if(exists("F_adjust2b")){
          F_adjust2b <- (F_adjust2b + global_adjuster)/2
        }else{
          F_adjust2b <- global_adjuster
        }
        F_adjust2a <- F_adjust2a*(Last_Mult2a-1)/(Last_Mult2a-median(Fmult2[adjusted_Rebuild_F_Rebuild]))
        F_adjust2b <- F_adjust2b*(Last_Mult2b-1)/(Last_Mult2b-median(Fmult2[adjusted_OFL_F_Rebuild]))
        if(!is.na(F_adjust2a)){
          if(is.infinite(F_adjust2a)){
            F_adjust2a <- global_adjuster
          }else if(F_adjust2a <= 0){
            F_adjust2a <- 0.5
          }
          if(is.infinite(F_adjust2a)){
            F_adjust2a <- global_adjuster
          }else if(F_adjust2a > 0){
            Fmult2a <- Fmult2[adjusted_Rebuild_F_Rebuild]
            if(length(Fmult2a[Fmult2a<1])>0){
              Fmult2a[Fmult2a<1] <- exp(log(Fmult2a[Fmult2a<1])*F_adjust2a)
            }
            if(length(Fmult2a[Fmult2a>1])>0){
              Fmult2a[Fmult2a>1] <- ((Fmult2a[Fmult2a>1]-1)*F_adjust2a+1)
            }
            Fmult2[adjusted_Rebuild_F_Rebuild] <- Fmult2a
        }else{
          F_adjust2a <- global_adjuster
        }}else{
          F_adjust2a <- global_adjuster
        }
        if(!is.na(F_adjust2b)){
          if(is.infinite(F_adjust2b)){
            F_adjust2b <- global_adjuster
          }else if(F_adjust2b <= 0){
            F_adjust2b <- 0.5
          }
          if(is.infinite(F_adjust2b)){
            F_adjust2b <- global_adjuster
          }else if(F_adjust2b > 0){
          Fmult2b <- Fmult2[adjusted_OFL_F_Rebuild]
          if(length(Fmult2b[Fmult2b<1])>0){
            Fmult2b[Fmult2b<1] <- exp(log(Fmult2b[Fmult2b<1])*F_adjust2b)
          }
          if(length(Fmult2b[Fmult2b>1])>0){
            Fmult2b[Fmult2b>1] <- ((Fmult2b[Fmult2b>1]-1)*F_adjust2b+1)
          }
          Fmult2[adjusted_OFL_F_Rebuild] <- Fmult2b
        }else{
          F_adjust2b <- global_adjuster
        }}else{
          F_adjust2b <- global_adjuster
        }
      }
    }
    Fmult1 <- rep(DepletionScale,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
	#Here the achieved catch fractions by fishing sector and year are calculated and compared relative 
	#to the target allocations. An adjustment multiplier is then computed to adjust fleet Fs closer to a
	#value expected to achieve the target allocations.
    if(FScale > 0){               
      if(n_groups>0){
        Catch_temp <- TimeFit3[,Catch_cols3]
        Catch_tot <- apply(Catch_temp[,unlist(fleets_by_group),drop=FALSE],1,sum)
        for(i in 1:n_groups){
          sort.mat <- matrix(NA, nrow = forecast[["Nforecastyrs"]]*length(seasons)*length(which(groups==i)), ncol = 2)
          sort.mat[,1] <- rep(1:forecast[["Nforecastyrs"]],length(seasons)*length(which(groups==i)))
          sort.mat[,2] <- rep(apply(Catch_temp[,fleets_by_group[[i]],drop=FALSE],1,sum)/Catch_tot,length(seasons)*length(which(groups==i)))
          sort.mat <- sort.mat[order(sort.mat[,1]),]
          Allocations[Allocations[,4]==i,6] <- sort.mat[,2]
        }
      }
      Fmult3 <- (0.5*(Allocations[,5]/Allocations[,6]-1)+1)
    }else{
      Fmult3 <- rep(1,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
    }
	  
    if(fitting_Rebuild==TRUE){
      if(Rebuild.Scale==0){
          Fmult3[adjusted_Rebuild_F_Rebuild] <- 1
        }
    }
    
    if(fitting_Fixed_Catch==FALSE){
      Fmult4 <- rep(1,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols))
    }
    
    
   
	#Adjust any multipliers of fixed catch values to 1 so that the 
    #search algorithm will consider them to have achieved their target	
    if(!is.null(fixed_ref)){
      Fmult1[fixed_ref] <- 1
      Fmult2[fixed_ref] <- 1
      Fmult3[fixed_ref] <- 1
      Fmult4[fixed_ref] <- 1
    }
    Comb_Mult <- Fmult1*Fmult2*Fmult3*Fmult4
    Comb_Mult[which(forecast_F[,4]>=max_F_limit & Comb_Mult>1)] <- 1
    
    #Record the previous adjustment values so they can be used to optimize 
    #step sizes to speed up target convergence	
    Last_Mult1 <- DepletionScale
    Last_Mult2 <- Fmult2
    if(!is.null(rebuild_ref)){
    Last_Mult2a <- median(Fmult2[adjusted_Rebuild_F_Rebuild])
    Last_Mult2b <- median(Fmult2[adjusted_OFL_F_Rebuild])
    }else{
      Last_Mult2a <- 1#rep(1,length(adjusted_F_OFL))
      Last_Mult2b <- 1#rep(1,length(adjusted_F_OFL))
    }
    Last_max_mult <- Curr_max_mult
    Curr_max_mult <- max(c(abs(1-Fmult1_raw),abs(1-Fmult2_raw),abs(1-Fmult3),abs(1-Fmult4)))
    Min_max_mult <- min(Min_max_mult,Curr_max_mult)
    
    
	#Plot out progess in achieving targets. This is primarily for diagnosis of a 
	#run that is failing to converge on an answer in a reasonable period of time.
    if(Make_plots==TRUE){
      try(
        {
          par(mar=c(4,3,3,2))
          col_options <- c("black","dark red","dark green","dark blue","orange","purple","red","green","blue","brown","pink","yellow",colors())
          point_options <- c(16,15,17,18,8,9,10,11,12,13,0,1,2,3,4,5,6,14,21,22,23,24,25,19,20)
          plot(Fmult1,xlab="year/season/fleet",ylab="Depletion Adjustment",col=rep(col_options[seq_along(F_cols)],forecast[["Nforecastyrs"]]*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),forecast[["Nforecastyrs"]]),main = paste0(method," loop = ",loop))
          plot(rep(F_adjust1,forecast[["Nforecastyrs"]]*length(seasons)*length(F_cols)),xlab="year/season/fleet",ylab="Depletion Optimization Adjustment",col=rep(col_options[seq_along(F_cols)],forecast[["Nforecastyrs"]]*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),forecast[["Nforecastyrs"]]),main = paste0(method," loop = ",loop))
          plot(Fmult2,xlab="year/season/fleet",ylab="F Adjustment",col=rep(col_options[seq_along(F_cols)],forecast[["Nforecastyrs"]]*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),forecast[["Nforecastyrs"]]),main = paste0(method," loop = ",loop))
          plot(F_adjust2,xlab="year/season/fleet",ylab="F Optimization Adjustment",col="red",pch=16,main = paste0(method," loop = ",loop))
          plot(Fmult3,xlab="year/season/fleet",ylab="Allocation Adjustment",col=rep(col_options[seq_along(F_cols)],forecast[["Nforecastyrs"]]*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),forecast[["Nforecastyrs"]]),main = paste0(method," loop = ",loop))
          plot(F_adjust2_2,xlab="year/season/fleet",ylab="Allocation Optimization Adjustment",col=rep(col_options[seq_along(F_cols)],forecast[["Nforecastyrs"]]*length(seasons)),pch=rep(sort(rep(point_options[seq_along(seasons)],length(F_cols))),forecast[["Nforecastyrs"]]),main = paste0(method," loop = ",loop))
        }
      )
    }
	#Check if all targets have been achieved and if so stop fitting
    if(max(abs(1-Fmult1_raw))>Depletion.Threshold | max(abs(1-Fmult2_raw))>Annual.F.Threshold | max(abs(1-Fmult3))>Allocation.Threshold | max(abs(1-Fmult4))>Annual.F.Threshold | abs(search_step)>Step.Threshold | loop < 2){keepFitting<-TRUE}else{keepFitting<-FALSE}
    if(FScale==0 & loop>2 & fitting_Fixed_Catch==FALSE){keepFitting<-FALSE}
  
    if(is.element(loop,c(1:30,seq(35,1000,5))) | global_adjuster<1){
      if(Messages == TRUE){
        message(paste0("Current loop = ",loop," ; still optimizing = ",keepFitting))
        message(paste0("Depletion optimization scaler = ",round(max(abs(1-Fmult1_raw)),6)," ; threshold <= ",Depletion.Threshold))
        message(paste0("Annual F optimization max scaler = ",round(max(abs(1-Fmult2_raw)),6)," ; threshold <= ",Annual.F.Threshold))
        message(paste0("Allocation optimization max scaler = ",round(max(abs(1-Fmult3)),6)," ; threshold <= ",Allocation.Threshold))
        message(paste0("Fixed catch optimization max scaler = ",round(max(abs(1-Fmult4)),6)," ; threshold <= ",Annual.F.Threshold))
        message(paste0("Step size optimization max scaler = ",search_step," ; threshold <= ",Step.Threshold))
        message(paste0("Global multiplier adjuster = ",global_adjuster))
        message(paste0("Best max multiplier = ",Min_max_mult))
      }
    }
    if(is.element(loop,c(50,100,200,500,1000,seq(1500,10000,500)))){
      Depletion.Threshold <- Depletion.Threshold*10
      Annual.F.Threshold <- Annual.F.Threshold*10
      Allocation.Threshold <- Allocation.Threshold*10
    }
	#Here we check that no Fs have been reduced to zero that need some catch
	#If that has occured repace the zero F with a small starting value 0.05 so that the 
	#search algorithm can act on it to achieve the true target value.
	#This is needed if the ABC loop was used to perform a zero catch run and then 
	#rebuild run is performed starting from those zero values
    zero_Fs <- which(forecast_F[,4]==0)
    increase_Fs <- which(Comb_Mult>1)
    if(length(zero_Fs)>0 & length(increase_Fs)>0){
     mod_Fs <- zero_Fs[which(is.element(zero_Fs,increase_Fs))]
      if(length(mod_Fs)>0){
        message("F values being reset from 0 to initial average this shouldn't happen so check results")
        forecast_F[mod_Fs,4] <- Forecast_catch_setup[mod_Fs,4]
        Comb_Mult[mod_Fs] <- runif(length(mod_Fs),0.95,1.05)
      }
    }
    forecast_F[,4] <- forecast_F[,4]*Comb_Mult
   } }
	#Now adjust the previous F values by the estimated multiplier to create a 
	#new estimate of the target Fs, make sure to overwrite any fixed catches 
	#with their original values.
    
    if(!is.null(fixed_ref)){
      forecast_F[fixed_ref,4] <- Fixed_catch_target[,4]
    }
    
    
    if(loop > 10){
      if(Curr_max_mult >= Last_max_mult){
        global_adjuster <- global_adjuster*0.8
        Min_max_mult <- 1.1*Min_max_mult
        Depletion.Threshold <- min(Depletion.Threshold*2,Min_max_mult)
        Annual.F.Threshold <- min(Annual.F.Threshold*2,Min_max_mult)
        Allocation.Threshold <- min(Allocation.Threshold*2,Min_max_mult)
      }
    }
    if(Curr_max_mult >= Last_max_mult){
      forecast_F <- last_forecast_F
      global_adjuster <- global_adjuster*0.95
      if(loop > 10){
        global_adjuster <- global_adjuster*0.85
        Min_max_mult <- 1.1*Min_max_mult
        Depletion.Threshold <- Depletion.Threshold*2
        Annual.F.Threshold <- Annual.F.Threshold*2
        Allocation.Threshold <- Allocation.Threshold*2
      }
    }else{
      last_forecast_F <- forecast[["ForeCatch"]]
    }
    Depletion.Threshold <- min(Depletion.Threshold,Min_max_mult)
    Annual.F.Threshold <- min(Annual.F.Threshold,Min_max_mult)
    Allocation.Threshold <- min(Allocation.Threshold,Min_max_mult)
    forecast[["ForeCatch"]] <- forecast_F
    #Write the modified forecast data out to a file and rerun projections
    unlink(paste0(getwd(),"/forecast.ss"))
    r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
    
    dir <- normalizePath(getwd())
    if(Run_in_MSE==FALSE){
      bin <- file.path(dir,SS_exe)
    }
    if (os == "unix") {
      system(
        paste0(
          "cd ", dir, ";", paste0(bin, " "),
          admb_options
        ),
        ignore.stdout = TRUE
      )
    } else {
      system(paste0(paste0(bin, " "), admb_options),
             invisible = TRUE, ignore.stdout = TRUE,
             show.output.on.console = FALSE
      )
    }
    Sys.sleep(0.05)
    
    #If all values have converged check if this is the OFL, ABC, or Rebuild loop
    if(keepFitting==FALSE){
      if(fitting_Fixed_Catch){
        if(Messages == TRUE){
          message(paste0("Constant Catch target ",CC_loop," of ",length(Const_Catch)," achieved for allocation loop ",allocation_loop," of ",N_allocation_scenarios))
        }
        if(Calc_Hessian==TRUE | Do_Pstar==TRUE){
          admb_options <- ""
          start$last_estimation_phase <- 10
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
          
          dir <- normalizePath(getwd())
          if(Run_in_MSE==FALSE){
            bin <- file.path(dir,SS_exe)
          }
          if (os == "unix") {
            system(
              paste0(
                "cd ", dir, ";", paste0(bin, " "),
                admb_options
              ),
              ignore.stdout = TRUE
            )
          } else {
            system(paste0(paste0(bin, " "), admb_options),
                   invisible = TRUE, ignore.stdout = TRUE,
                   show.output.on.console = FALSE
            )
          }
          Sys.sleep(0.05)
          
          resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
          
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("Param_CC_",CC_loop)]] <- resultsFit$parameters[,c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("F_sd_CC_",CC_loop)]] <- resultsFit$parameters[grep("F_fleet",resultsFit$parameters$Label),c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("SSB_sd_CC_",CC_loop)]] <- resultsFit$derived_quants[grep("SSB_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("Recr_sd_CC_",CC_loop)]] <- resultsFit$derived_quants[grep("Recr_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("ForeCatchF_CC_",CC_loop)]] <- resultsFit$derived_quants[grep("F_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("ForeCatchF_CC_",CC_loop)]] <- projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("ForeCatchF_CC_",CC_loop)]][-grep("annF_",projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("ForeCatchF_CC_",CC_loop)]][,"Label"]),]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("ForeCatchDead_CC_",CC_loop)]] <- resultsFit$derived_quants[grep("ForeCatch_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][[paste0("ForeCatchRetained_CC_",CC_loop)]] <- resultsFit$derived_quants[grep("ForeCatchret_",resultsFit$derived_quants$Label),1:3]
          
          start$last_estimation_phase <- 0
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
        }
        admb_options <- "-nohess"
        if(allocation_loop==1){
          if(dir.exists(paste0(Assessment_dir,"/CC_",CC_loop))){
            unlink(paste0(Assessment_dir,"/CC_",CC_loop),recursive = TRUE)
          }
          dir.create(paste0(Assessment_dir,"/CC_",CC_loop))
          temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
          file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/CC_",CC_loop,"/",temp.files))
          
          CC_ForeCatch_1 <- forecast_F
        }
        dir.create(paste0(Assessment_dir,"/CC_",CC_loop,"/Allocation_run_",allocation_loop))
        temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
        file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/CC_",CC_loop,"/Allocation_run_",allocation_loop,"/",temp.files))
        
        if(allocation_loop<N_allocation_scenarios){
          allocation_loop <- allocation_loop+1
          keepFitting <- TRUE
        }else if(length(Const_Catch)>CC_loop){
          allocation_loop <- 1
          CC_loop <- CC_loop + 1
          Catch_Target <- rep(Const_Catch[CC_loop],forecast[["Nforecastyrs"]])
          if(Catch_trunc > 0){
            Catch_Target[(forecast$Nforecastyrs-c((Catch_trunc-1):0))] <- 0
          }
          keepFitting <- TRUE
          
          #Return ForeCatch F's to the first allocation estimate from last 
          #constant catch this should be closer to the real solution.
          
          forecast_F <- CC_ForeCatch_1 
          forecast[["ForeCatch"]] <- forecast_F
          unlink(paste0(getwd(),"/forecast.ss"))
          r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
        }
        Allocations <- Allocation_tracker[[allocation_loop]]  
      }else if(fitting_Benchmark==TRUE){
        if(Messages == TRUE){
          message(paste0("Benchmark target achieved for allocation loop ",allocation_loop," of ",N_allocation_scenarios))
        }
        if(Calc_Hessian==TRUE | Do_Pstar==TRUE){
          admb_options <- ""
          start$last_estimation_phase <- 10
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
          
          dir <- normalizePath(getwd())
          if(Run_in_MSE==FALSE){
            bin <- file.path(dir,SS_exe)
          }
          if (os == "unix") {
            system(
              paste0(
                "cd ", dir, ";", paste0(bin, " "),
                admb_options
              ),
              ignore.stdout = TRUE
            )
          } else {
            system(paste0(paste0(bin, " "), admb_options),
                   invisible = TRUE, ignore.stdout = TRUE,
                   show.output.on.console = FALSE
            )
          }
          Sys.sleep(0.05)
          
          
          resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
          
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Param_Benchmark"]] <- resultsFit$parameters[,c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_sd_Benchmark"]] <- resultsFit$parameters[grep("F_fleet",resultsFit$parameters$Label),c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_sd_Benchmark"]] <- resultsFit$derived_quants[grep("SSB_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recr_sd_Benchmark"]] <- resultsFit$derived_quants[grep("Recr_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Benchmark"]] <- resultsFit$derived_quants[grep("F_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Benchmark"]] <- projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Benchmark"]][-grep("annF_",projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Benchmark"]][,"Label"]),]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchDead_Benchmark"]] <- resultsFit$derived_quants[grep("ForeCatch_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchRetained_Benchmark"]] <- resultsFit$derived_quants[grep("ForeCatchret_",resultsFit$derived_quants$Label),1:3]
          
          start$last_estimation_phase <- 0
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
        }
        admb_options <- "-nohess"
        #Write out the Benchmark results to a new folder (replace any old folder that exists)
        if(allocation_loop==1){
          if(dir.exists(paste0(Assessment_dir,"/Benchmark_target"))){
            unlink(paste0(Assessment_dir,"/Benchmark_target"),recursive = TRUE)
          }
          dir.create(paste0(Assessment_dir,"/Benchmark_target"))
          temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
          file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/Benchmark_target/",temp.files))
          
          Benchmark_ForeCatch_1 <- forecast_F
        }
        dir.create(paste0(Assessment_dir,"/Benchmark_target","/Allocation_run_",allocation_loop))
        temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
        file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/Benchmark_target/","/Allocation_run_",allocation_loop,"/",temp.files))
        
        if(allocation_loop<N_allocation_scenarios){
          allocation_loop <- allocation_loop+1
          keepFitting <- TRUE
        }else {
          if(Messages == TRUE){
            message(paste0("Beginning OFL optimization"))
          }
          allocation_loop <- 1
          fitting_Benchmark <- FALSE
          fitting_OFL <- TRUE
          par(mfrow=c(4,2))
          loop <- 0
          method <- "OFL"
          keepFitting <- TRUE
          forecast$fcast_rec_option <- OFL_recruit_setting 
          forecast_F <- Benchmark_ForeCatch_1 
          forecast[["ForeCatch"]] <- forecast_F
          r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
        }
        Allocations <- Allocation_tracker[[allocation_loop]]  
      }else if(fitting_OFL==TRUE){
        if(Messages == TRUE){
          message(paste0("OFL target achieved for allocation loop ",allocation_loop," of ",N_allocation_scenarios))
        }
        if(Calc_Hessian==TRUE){
          admb_options <- ""
          start$last_estimation_phase <- 10
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
          
          dir <- normalizePath(getwd())
          if(Run_in_MSE==FALSE){
            bin <- file.path(dir,SS_exe)
          }
          if (os == "unix") {
            system(
              paste0(
                "cd ", dir, ";", paste0(bin, " "),
                admb_options
              ),
              ignore.stdout = TRUE
            )
          } else {
            system(paste0(paste0(bin, " "), admb_options),
                   invisible = TRUE, ignore.stdout = TRUE,
                   show.output.on.console = FALSE
            )
          }
          Sys.sleep(0.05)
          
          
          resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
          
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Param_OFL"]] <- resultsFit$parameters[,c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_sd_OFL"]] <- resultsFit$parameters[grep("F_fleet",resultsFit$parameters$Label),c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_sd_OFL"]] <- resultsFit$derived_quants[grep("SSB_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recr_sd_OFL"]] <- resultsFit$derived_quants[grep("Recr_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_OFL"]] <- resultsFit$derived_quants[grep("F_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_OFL"]] <- projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_OFL"]][-grep("annF_",projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_OFL"]][,"Label"]),]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchDead_OFL"]] <- resultsFit$derived_quants[grep("ForeCatch_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchRetained_OFL"]] <- resultsFit$derived_quants[grep("ForeCatchret_",resultsFit$derived_quants$Label),1:3]
          
          start$last_estimation_phase <- 0
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
        }
        admb_options <- "-nohess"
        
        if(allocation_loop==1){
          #Write out the OFL results to a new folder (replace any old folder that exists)
          if(dir.exists(paste0(Assessment_dir,"/OFL_target"))){
            unlink(paste0(Assessment_dir,"/OFL_target"),recursive = TRUE)
          }
          dir.create(paste0(Assessment_dir,"/OFL_target"))
          temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
          file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/OFL_target/",temp.files))
          
          OFL_ForeCatch_1 <- forecast_F
        }
        dir.create(paste0(Assessment_dir,"/OFL_target","/Allocation_run_",allocation_loop))
        temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
        file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/OFL_target/","/Allocation_run_",allocation_loop,"/",temp.files))
        
        if(allocation_loop<N_allocation_scenarios){
          allocation_loop <- allocation_loop+1
          keepFitting <- TRUE
        }else {
          fitting_OFL <- FALSE
          par(mfrow=c(4,2))
          forecast_F <- OFL_ForeCatch_1 
          forecast[["ForeCatch"]] <- forecast_F
          allocation_loop <- 1
          if(!is.null(ABC_Fraction)){
            if(Messages == TRUE){
              message(paste0("Beginning ABC optimization"))
            }
            #If in the OLF loop then reset to keep fitting and change the target from OFL to ABC
            loop <- 0
            method <- "ABC"
            forecast$fcast_rec_option <- ABC_recruit_setting
            r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
            keepFitting <- TRUE
            fitting_ABC <- TRUE
          }else if(Calc_F0==TRUE){
            if(Messages == TRUE){
              message(paste0("Beginning F0 optimization"))
            }
            loop <- 0
            method <- "F0"
            forecast$fcast_rec_option <- F0_recruit_setting
            r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
            keepFitting <- TRUE
            fitting_F0 <- TRUE
          }else if(!is.null(Rebuild_yr)){
            if(Messages == TRUE){
              message(paste0("Beginning Rebuild optimization"))
            }
            #If in the OLF loop then reset to keep fitting and change the target from OFL to Rebuild
            loop <- 0
            method <- "Rebuild"
            forecast$fcast_rec_option <- Rebuild_recruit_setting
            r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
            keepFitting <- TRUE
            fitting_Rebuild <- TRUE
          }
        }
        Allocations <- Allocation_tracker[[allocation_loop]]  
      }else if(fitting_ABC==TRUE){
        if(Messages == TRUE){
          message(paste0("ABC target achieved for allocation loop ",allocation_loop," of ",N_allocation_scenarios))
        }
        if(Calc_Hessian==TRUE){
          admb_options <- ""
          start$last_estimation_phase <- 10
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
          
          dir <- normalizePath(getwd())
          if(Run_in_MSE==FALSE){
            bin <- file.path(dir,SS_exe)
          }
          if (os == "unix") {
            system(
              paste0(
                "cd ", dir, ";", paste0(bin, " "),
                admb_options
              ),
              ignore.stdout = TRUE
            )
          } else {
            system(paste0(paste0(bin, " "), admb_options),
                   invisible = TRUE, ignore.stdout = TRUE,
                   show.output.on.console = FALSE
            )
          }
          Sys.sleep(0.05)
          
          
          resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
          
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Param_ABC"]] <- resultsFit$parameters[,c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_sd_ABC"]] <- resultsFit$parameters[grep("F_fleet",resultsFit$parameters$Label),c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_sd_ABC"]] <- resultsFit$derived_quants[grep("SSB_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recr_sd_ABC"]] <- resultsFit$derived_quants[grep("Recr_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_ABC"]] <- resultsFit$derived_quants[grep("F_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_ABC"]] <- projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_ABC"]][-grep("annF_",projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_ABC"]][,"Label"]),]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchDead_ABC"]] <- resultsFit$derived_quants[grep("ForeCatch_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchRetained_ABC"]] <- resultsFit$derived_quants[grep("ForeCatchret_",resultsFit$derived_quants$Label),1:3]
          
          start$last_estimation_phase <- 0
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
        }
        admb_options <- "-nohess"
        
        if(allocation_loop==1){
          #Write out the OFL results to a new folder (replace any old folder that exists)
          if(dir.exists(paste0(Assessment_dir,"/ABC_target"))){
            unlink(paste0(Assessment_dir,"/ABC_target"),recursive = TRUE)
          }
          dir.create(paste0(Assessment_dir,"/ABC_target"))
          temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
          file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/ABC_target/",temp.files))
          
          ABC_ForeCatch_1 <- forecast_F
        }
        dir.create(paste0(Assessment_dir,"/ABC_target","/Allocation_run_",allocation_loop))
        temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
        file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/ABC_target/","/Allocation_run_",allocation_loop,"/",temp.files))
        
        if(allocation_loop<N_allocation_scenarios){
          allocation_loop <- allocation_loop+1
          keepFitting <- TRUE
        }else {
          fitting_ABC <- FALSE
          par(mfrow=c(4,2))
          forecast_F <- ABC_ForeCatch_1 
          forecast[["ForeCatch"]] <- forecast_F
          allocation_loop <- 1
          if(Calc_F0==TRUE){
            if(Messages == TRUE){
              message(paste0("Beginning F0 optimization"))
            }
            loop <- 0
            method <- "F0"
            forecast$fcast_rec_option <- F0_recruit_setting
            r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
            keepFitting <- TRUE
            fitting_F0 <- TRUE
          }else if(!is.null(Rebuild_yr)){
            if(Messages == TRUE){
              message(paste0("Beginning Rebuild optimization"))
            }
            #If in the ABC loop then reset to keep fitting and change the target from ABC to Rebuild
            loop <- 0
            method <- "Rebuild"
            forecast$fcast_rec_option <- Rebuild_recruit_setting
            r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
            keepFitting <- TRUE
            fitting_Rebuild <- TRUE
          }
        }
        Allocations <- Allocation_tracker[[allocation_loop]]  
        
      }else if(fitting_F0==TRUE){
        if(Messages == TRUE){
          message(paste0("F0 target achieved for allocation loop ",allocation_loop," of ",N_allocation_scenarios))
        }
        if(Calc_Hessian==TRUE){
          admb_options <- ""
          start$last_estimation_phase <- 10
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
          
          dir <- normalizePath(getwd())
          if(Run_in_MSE==FALSE){
            bin <- file.path(dir,SS_exe)
          }
          if (os == "unix") {
            system(
              paste0(
                "cd ", dir, ";", paste0(bin, " "),
                admb_options
              ),
              ignore.stdout = TRUE
            )
          } else {
            system(paste0(paste0(bin, " "), admb_options),
                   invisible = TRUE, ignore.stdout = TRUE,
                   show.output.on.console = FALSE
            )
          }
          Sys.sleep(0.05)
          
          
          resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
          
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Param_F0"]] <- resultsFit$parameters[,c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_sd_F0"]] <- resultsFit$parameters[grep("F_fleet",resultsFit$parameters$Label),c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_sd_F0"]] <- resultsFit$derived_quants[grep("SSB_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recr_sd_F0"]] <- resultsFit$derived_quants[grep("Recr_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_F0"]] <- resultsFit$derived_quants[grep("F_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_F0"]] <- projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_F0"]][-grep("annF_",projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_F0"]][,"Label"]),]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchDead_F0"]] <- resultsFit$derived_quants[grep("ForeCatch_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchRetained_F0"]] <- resultsFit$derived_quants[grep("ForeCatchret_",resultsFit$derived_quants$Label),1:3]
          
          start$last_estimation_phase <- 0
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
        }
        admb_options <- "-nohess"
        
        if(allocation_loop==1){
          #Write out the OFL results to a new folder (replace any old folder that exists)
          if(dir.exists(paste0(Assessment_dir,"/F0_target"))){
            unlink(paste0(Assessment_dir,"/F0_target"),recursive = TRUE)
          }
          dir.create(paste0(Assessment_dir,"/F0_target"))
          temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
          file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/F0_target/",temp.files))
          
          F0_ForeCatch_1 <- forecast_F
        }
        dir.create(paste0(Assessment_dir,"/F0_target","/Allocation_run_",allocation_loop))
        temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
        file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/F0_target/","/Allocation_run_",allocation_loop,"/",temp.files))
        
        if(allocation_loop<N_allocation_scenarios){
          allocation_loop <- allocation_loop+1
          keepFitting <- TRUE
        }else {
          fitting_F0 <- FALSE
          par(mfrow=c(4,2))
          forecast_F <- OFL_ForeCatch_1 
          forecast[["ForeCatch"]] <- forecast_F
          allocation_loop <- 1
          if(!is.null(Rebuild_yr)){
            if(Messages == TRUE){
              message(paste0("Beginning Rebuild optimization"))
            }
            loop <- 0
            method <- "Rebuild"
            forecast$fcast_rec_option <- Rebuild_recruit_setting
            r4ss::SS_writeforecast(mylist=forecast,overwrite = TRUE, verbose = Verbose)
            keepFitting <- TRUE
            fitting_Rebuild <- TRUE
          }
        }
        
        Allocations <- Allocation_tracker[[allocation_loop]]  
      }else if(fitting_Rebuild==TRUE){
        if(Messages == TRUE){
          message(paste0("Rebuild target achieved for allocation loop ",allocation_loop," of ",N_allocation_scenarios))
        }
        if(Calc_Hessian==TRUE){
          admb_options <- ""
          start$last_estimation_phase <- 10
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
          
          dir <- normalizePath(getwd())
          if(Run_in_MSE==FALSE){
            bin <- file.path(dir,SS_exe)
          }
          if (os == "unix") {
            system(
              paste0(
                "cd ", dir, ";", paste0(bin, " "),
                admb_options
              ),
              ignore.stdout = TRUE
            )
          } else {
            system(paste0(paste0(bin, " "), admb_options),
                   invisible = TRUE, ignore.stdout = TRUE,
                   show.output.on.console = FALSE
            )
          }
          Sys.sleep(0.05)
          
          
          resultsFit <- r4ss::SS_output(dir=getwd(),covar=FALSE, verbose = Verbose, printstats= Verbose)
          
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Param_Rebuild"]] <- resultsFit$parameters[,c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["F_sd_Rebuild"]] <- resultsFit$parameters[grep("F_fleet",resultsFit$parameters$Label),c(2,3,11,12)]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["SSB_sd_Rebuild"]] <- resultsFit$derived_quants[grep("SSB_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["Recr_sd_Rebuild"]] <- resultsFit$derived_quants[grep("Recr_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Rebuild"]] <- resultsFit$derived_quants[grep("F_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Rebuild"]] <- projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Rebuild"]][-grep("annF_",projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchF_Rebuild"]][,"Label"]),]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchDead_Rebuild"]] <- resultsFit$derived_quants[grep("ForeCatch_",resultsFit$derived_quants$Label),1:3]
          projection_results[[paste0("Allocation_run_",allocation_loop)]][["ForeCatchRetained_Rebuild"]] <- resultsFit$derived_quants[grep("ForeCatchret_",resultsFit$derived_quants$Label),1:3]
          
          start$last_estimation_phase <- 0
          r4ss::SS_writestarter(mylist=start,overwrite = TRUE, verbose = Verbose)
        }
        admb_options <- "-nohess"
        
        if(allocation_loop==1){
          #Write out the OFL results to a new folder (replace any old folder that exists)
          if(dir.exists(paste0(Assessment_dir,"/Rebuild_target"))){
            unlink(paste0(Assessment_dir,"/Rebuild_target"),recursive = TRUE)
          }
          dir.create(paste0(Assessment_dir,"/Rebuild_target"))
          temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
          file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/Rebuild_target/",temp.files))
          
          Rebuild_ForeCatch_1 <- forecast_F
        }
        dir.create(paste0(Assessment_dir,"/Rebuild_target","/Allocation_run_",allocation_loop))
        temp.files <- list.files(path=paste0(Assessment_dir,"/Working_dir"))
        file.copy(from=paste0(Assessment_dir,'/Working_dir/',temp.files),to=paste0(Assessment_dir,"/Rebuild_target/","/Allocation_run_",allocation_loop,"/",temp.files))
        
        if(allocation_loop<N_allocation_scenarios){
          allocation_loop <- allocation_loop+1
          keepFitting <- TRUE
        }else {
          fitting_Rebuild <- FALSE
          par(mfrow=c(4,2))
          forecast_F <- Rebuild_ForeCatch_1 
          forecast[["ForeCatch"]] <- forecast_F
          allocation_loop <- 1
        }
        Allocations <- Allocation_tracker[[allocation_loop]]  
        
      }
    }
  }
  setwd(oldwd)
  if(Messages == TRUE){
    message(paste0("Projection run complete"))
  }
  if(return_warnings==TRUE){
    warning(paste0(combinded_warnings))
  }
  
  invisible(projection_results)
}
