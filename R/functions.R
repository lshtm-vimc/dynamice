# functions.R
# Functions for DynaMICE (Dynamic Measles Immunisation Calculation Engine)

# Measles vaccine impact model
#   To estimate the health impact of measles vaccination for a given set of 
#   countries.

# ------------------------------------------------------------------------------
# functions
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# create_life_expectancy_remaining_full
#
# create remaining life expectancy file for each year across all age intervals
# ------------------------------------------------------------------------------
#   uses folder - input  -- "input/demographicdata2019/"
#   uses folder - output -- "input/demographicdata2019/"
# ------------------------------------------------------------------------------
create_life_expectancy_remaining_full <- function () {
  
  # read data file for remaining life expectancy
  lexp_remain       <- fread ("input/demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both.csv")
  
  lexp_remain_full <- copy (lexp_remain)
  
  # Since year is 5-year intervals, add data for all in-between years
  for (i in 1:4) {
    
    dt <- copy (lexp_remain) 
    dt [, year := year + i]
    
    # add table to full table of remaining life expectancy 
    lexp_remain_full <- rbindlist (list (lexp_remain_full, dt), 
                                   use.names = T)
  }
  
  # add data for 2100 for remaining life expectancy
  # copy values for year 2095 to year 2100 
  lexp_remain_2100 <- lexp_remain [year == 2095]
  lexp_remain_2100 [, year := 2100]
  
  lexp_remain_full <- rbindlist (list (lexp_remain_full, lexp_remain_2100), 
                                 use.names = T)
  
  # save file
  fwrite (lexp_remain_full, 
          file = "input/demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both_full.csv")
  
} # end of function -- create_life_expectancy_remaining_full
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create_vaccine_coverage_routine_sia
#
# create 2 vaccine coverage files per scenario for routine and 
# SIA (supplementary immunisation activities) from VIMC vaccine coverage files
# ------------------------------------------------------------------------------
create_vaccine_coverage_routine_sia <- function (vaccine_coverage_folder    = "", 
                                                 coverage_prefix            = "",
                                                 touchstone                 = "",
                                                 antigen                    = "",
                                                 scenario_name,
                                                 vaccine_coverage_subfolder = "") {
  
  # vaccine coverage file
  vaccine_coverage_file <- paste0 (vaccine_coverage_folder, 
                                   coverage_prefix, 
                                   touchstone, 
                                   antigen,
                                   scenario_name, 
                                   ".csv")
  
  # separate vaccine coverage file for routine immunisation
  routine_coverage_file <- paste0 (vaccine_coverage_folder, 
                                   vaccine_coverage_subfolder,
                                   "routine_",
                                   scenario_name, 
                                   ".csv")
  
  # separate vaccine coverage file for SIA (supplementary immunisation activities)
  sia_coverage_file <- paste0 (vaccine_coverage_folder, 
                               vaccine_coverage_subfolder,
                               "sia_",
                               scenario_name, 
                               ".csv")
  
  # read vaccine coverage data file
  sia <- fread (file             = vaccine_coverage_file, 
                stringsAsFactors = F, 
                na.strings       = "<NA>")
  
  # extract routine vaccination coverage
  routine <- sia [activity_type != "campaign",]
  
  # only select campaigns
  sia2 <- sia [activity_type == "campaign" & ( (!is.na(target) & target != 0) | (!is.na(coverage) ) ), ]
  
  # ----------------------------------------------------------------------------
  # set start age to fraction of a year for campaigns starting at ages 
  # less than 1 year, that is, 6 months, 9 months, etc will become 0.5, 0.75, etc
  # If age_first = 1 and age_range_verbatim = default or other text, then set start age to 0.5 (6 months).

  sia2 [                 , age_first := as.double (age_first)]
  # sia2 [age_first == 1, age_first := as.double (str_extract (age_range_verbatim, "\\d+")) / 12]
  # age_range_verbatim "1-14 Y" should not interpreted as 1 month
  sia2 [age_first == 1 & (str_extract (age_range_verbatim, "\\d+") != "1"), age_first := as.double (str_extract (age_range_verbatim, "\\d+")) / 12]
  sia2 [is.na (age_first), age_first := 0.5]
  
  sia2 [age_first == 0,    age_first := as.double (str_extract (age_range_verbatim, "\\d+")) / 12]
  sia2 [is.na (age_first), age_first := 0.5]
  # ----------------------------------------------------------------------------

  
  # remove unused columns
  discard.cols <- c ("scenario", "set_name", "gavi_support", "activity_type")
  sia2 <- sia2 [, .SD, .SDcols = -discard.cols]
  
  # remove all whitespace from age_range_verbatim, put all in lowercase
  sia2$age_range_verbatim <- tolower (gsub("[[:space:]]", "", sia2$age_range_verbatim))
  
  # use more inprecise values if age range verbatim is empty
  temp.cols2 <- c("t_combined", "t_l_a0", "t_l_a1", "t_n_a0", "t_n_a1")
  
  # age_first and age_last is what is used
  sia2 [, (temp.cols2) := list (0, "y", "y", age_first, age_last)]
  
  # calculate values for a0 and a1
  # for age <= 3 years, age is in weeks
  # for age >= 3 years, (156 weeks for upto age 3) + (remaining years above 3)
  #      age 1 year ~ 52; 2 year ~ 104; 3 year ~ 156; 4 year ~ 157; 5 year ~ 158
  find.a <- function (x) {
    t0 <- ifelse (x <= 3, round (x * 52), round ((x - 3) + (3 * 52) ) )
    return (t0)
  }
  
  # set start and end ages for SIA
  sia2 [, `:=` (a0 = find.a(t_n_a0), a1 = find.a(t_n_a1))]
  
  # age 0 ~ 1 week
  sia2 [a0 == 0, a0 := 1]  
  sia2 [a1 == 0, a1 := 1]
  
  # not sure what next 2 lines are for? (check)
  sia2 [t_l_a0 != "y",]
  sia2 [t_l_a1 != "y",]
  
  # remove temporary columns
  sia2 <- sia2 [, .SD, .SDcols = -temp.cols2]
  
  # ----------------------------------------------------------------------------
  # FROM GUIDELINES FOR 2019 RUNS: 
  # *coverage* shows the level of vaccination coverage, usually ranging from 0 (0%)
  # to 1 (100%). In some cases, this value may be greater than 1. For example, if a
  # campaign originally targets 1 million people but ends up vaccinating 1.1
  # million people, the coverage would be shown as 1.1 (equating to 110%). 
  #
  # Coverage and target population are now always specified at a national level. 
  # 
  # For example, where a campaign targets all ages in Region A (population 1,000,000)
  # and achieves 90% coverage, and where the population of the whole country is
  # 5,000,000, the coverage would appear on Montagu as 0.18 (18%) and the target
  # population as 5,000,000. 
  # 
  # This way of specifying coverage and target population
  # applies in both past and future years. *target* is the number of individuals in
  # the target population. 
  #
  # This is always shown for campaigns, and is now specified
  # at a national level. (See ‘coverage’ section above.) 
  #
  # For routine, target is shown as NA, which means you should assume the target 
  # population matches the population shown in the demographic data downloads 
  # for the corresponding ages (age_first and age_last). 
  # ----------------------------------------------------------------------------
  
  # infer reached from coverage
  sia2 [, reached := as.numeric(target) * as.numeric(coverage)]
  
  
  # add additional columns to make file similar to original
  sia2 [, c("Geography", "Gavi73", "Gavi-supported", "Activity","Extent")] <- NA
  
  # reorder columns to make file similar to original
  sia2 <- sia2 [ , c("Geography", "country_code", "vaccine", "Gavi73", "Gavi-supported",
                     "Activity", "vaccine", "year", "a0", "a1", 
                     "Extent", "target", "reached", "coverage") ]
  
  sia2 [, vaccine := NULL]
  
  # write vaccine coverage data for routine vaccination and SIAs
  fwrite (x    = routine, 
          file = routine_coverage_file)
  
  fwrite (x    = sia2, 
          file = sia_coverage_file)
  
} # end of function -- create_vaccine_coverage_routine_sia
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# writelog
#
# small function to write message a to logfile logname
# ------------------------------------------------------------------------------
writelog <- function (logname, 
                      x) {
  write (
    paste0 (
      format (Sys.time(), "%Y/%m/%d %H:%M:%S"),
      " ",
      x
    ),
    file   = logname,
    append = TRUE
  )
  
} # end of function -- writelog
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# logwarning
#
# prints a message, and writes to the logfile
# ------------------------------------------------------------------------------
logwarning <- function (logname, 
                        x) {
  
  writelog (logname, x)
  print (x)
  
} # end of function -- logwarning
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# expandMatrix
#
# easily expand a Matrix (i.e. duplicate columns and rows in an AxB matrix - and rescale where applicable)
# easy and transparant way to redistribute contact matrix - set default rescale to FALSE
# ------------------------------------------------------------------------------
expandMatrix <- function (A, 
                          expand_rows  = 1, 
                          expand_cols  = 1, 
                          rescale_rows = F, 
                          rescale_cols = F) {
  
  if(!is.matrix(A)){
    stop("A is not a matrix")
  }
  
  matvals <- numeric(0)
  rows <- nrow(A)
  cols <- ncol(A)
  
  for(c in 1:cols) {
    matvals <- c(
      matvals,
      rep(
        A[,c],
        expand_cols,
        each = expand_rows
      )
    )
  }
  
  B <- matrix (matvals, 
               nrow = rows * expand_rows, 
               ncol = cols * expand_cols)
  
  if(rescale_rows & rescale_cols){
    B <- B/(expand_rows*expand_cols) 
  } else if(rescale_rows){
    B <- B/expand_rows
  } else if(rescale_cols){
    B <- B/expand_cols
  }
  
  return (B)
  
} # end of function -- expandMatrix
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# updateProgress
#
# update file showing progress (temporarily disabled)
# ------------------------------------------------------------------------------
updateProgress <- function (country,
                            cnt,
                            runs,
                            r,
                            status) {
  ## wait until gavi_progress is yours (so cannot be overwritten by two parallel scripts at the same time)
  # while(
  #	 !(file.rename("gavi_progress",paste0("gavi_progress_locked_",cnt,"_",r)))
  # ){}
  ## update line
  # progress <- readLines("gavi_progress_locked",-1)
  # progress[((cnt-1)*runs+r)] <- paste0(
  #	 country,
  #	 " ",
  #	 paste0(c(rep(0,(nchar(runs)-nchar(r))),r),collapse=""),
  #	 " ",
  #	 status
  # )
  # writeLines(progress,"gavi_progress_locked")
  # file.rename(paste0("gavi_progress_locked_",cnt,"_",r),"gavi_progress")
  
} # end of function -- updateProgress
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# timelinessCov
#
# timeliness of vaccination
# ------------------------------------------------------------------------------
timelinessCov <- function (tim_shape,
                           target_age,
                           max_cov) {
  
  coverage_byage_cumulative <- c(
    rep(0,(target_age-1)),
    rep(max_cov,(254-target_age+1))
  )
  
  beta_function <- cumsum(
    dbeta(
      (1:52/52),
      4,
      tim_shape
    )
  )/max(
    cumsum(
      dbeta(
        (1:52/52),
        4,
        tim_shape
      )
    )
  )
  
  coverage_byage_cumulative[39:(39+51)] <- max_cov*beta_function
  coverage_byage=rep (0, 254)
  coverage_byage [2:length(coverage_byage_cumulative)] <- (
    coverage_byage_cumulative [2:length(coverage_byage_cumulative)]
    -coverage_byage_cumulative [1:(length(coverage_byage_cumulative)-1)]
  )/(
    1 - coverage_byage_cumulative [1:(length(coverage_byage_cumulative)-1)]
  )
  
  coverage_byage [which (is.na (coverage_byage))] <- 1
  
  return (coverage_byage)
  
} # end of function -- timelinessCov
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# runCountry
#
# Run the measles model for a given set of countries
# ------------------------------------------------------------------------------
# Changes 2019:
# 1) expanded matrix is rescaled to keep the total number of
# contacts in weekly age groups correct (checked - R0 of rescaled matrix is now
# the same as prior to expanding it to weekly ages)
# 2) projects the contact matrix to represent country's demography
# ------------------------------------------------------------------------------
runCountry <- function (
  #variables specific for loop
  ii,
  countries,
  years,
  
  #infection dynamic variables
  vaccination,
  using_sia,
  dinf,
  gamma,
  r0_basic,
  amplitude,
  take,
  degree,
  sia.method,
  
  #input data
  coverage_routine,
  coverage_sia,
  timeliness,
  contact,
  rnought,
  population,
  lexp,
  cfr,
  
  # dynaMICE model options
  tstep,
  save.scenario,
  foldername,
  measles_model,
  debug_model,
  debug_spinup,
  debug_age,
  debug_compartments,
  debug_country,
  debug_relative,
  debug_timepoints,
  
  # PSA options
  r,
  runs,
  psa,
  psa_var,  # used only when psa > 0
  
  #additional options
  process_results,
  run_model,
  remove_files = TRUE
) {
  
  iso3 <- countries[ii]
  
  # temporarily assign 3-letter-ISO code to Kosovo until Kosovo is assigned official ISO3-code
  if (iso3 == "XK"){
    fortran_country_code <- "XKX"
  } else {
    fortran_country_code <- iso3
  }
  
  if(psa > 0) {
    p <- r
    if(p<10){
      p <- paste0("00",p)
    } else if(p<100){
      p <- paste0("0",p)
    }
    psadir <- paste0("run",p,"/")
  } else {
    psadir <- ""
  }
  
  if (run_model) {
    
    if (using_sia == 1) {
      
      if (nrow (coverage_sia [country_code == iso3]) > 0) {
        
        sia_a0       <- coverage_sia [country_code == iso3, a0]
        sia_year     <- coverage_sia [country_code == iso3, year]
        sia_a1       <- coverage_sia [country_code == iso3, a1]
        sia_coverage <- coverage_sia [country_code == iso3, coverage]
        
        sia_coverage [which (sia_coverage>1)] <- 0.95
        
      } else {
        
        sia_a0 			 <- 0
        sia_a1 			 <- 0
        sia_coverage <- 0
        sia_year 		 <- years[1]
      }
    } else {
      
      sia_a0 			 <- 0
      sia_a1 			 <- 0
      sia_coverage <- 0
      sia_year 		 <- years[1]
    }
    
    # t_cov is x year of model (model is ran for t_cov timesteps)
    t_cov <- length(years) * tstep
    t_end <- t_cov * 2
    
    # t_run: for each year that is modelled,  get timepoint at which the year starts
    t_run <- t_end - (length(years):1) * tstep + 1
    
    # get timepoints for year sia, assuming that year 1 of interest is the year 2000
    t_sia <- (t_run[1] + 1000 * (sia_year - years[1])) + c(1:length(sia_year))
    
    if (using_sia == 0) {
      t_sia <- t_sia*-1
    }
    
    print ("Generating data for model...")
    writelog ("gavi_log", paste0 (iso3, "; Run ",r,"/",runs,"; Start generating data"))
    updateProgress (iso3, ii, runs, r, 1)
    
    # Vaccine parameters
    if (psa > 0) {
      take <- as.numeric (psa_var [r, c("take1_input",
                                        "take2_input",
                                        "take3_input")])
    }
    
    # vaccine take
    take1 <- c (take[1], 0)[sia.method]    # vaccine take (dose 1, before age 1)
    take2 <- c (take[2], 0)[sia.method]    # vaccine take (dose 1, after age 1)
    take3 <- c (take[3], 0)[sia.method]    # vaccine take (dose 2)
    
    # vaccine degree
    degree1 <- c(0, degree[1])[sia.method]  # vaccine degree (dose 1, before age 1)
    degree2 <- c(0, degree[2])[sia.method]  # vaccine degree (dose 1, after age 1)
    degree3 <- c(0, degree[3])[sia.method]  # vaccine degree (dose 2)
    
    # country_specific contact matrix
    r0_target     <- rnought [country_code == iso3, r0]
    q             <- r0_target / r0_basic  # proportionality factor (infectivity, underreporting)
    gamma         <- 1 / ( dinf * tstep / 365 )  # rate of losing infection
    contact_day   <- contact * q
    contact_tstep <- contact_day * (365 / tstep)
    r0_tstep      <- Re (eigen (contact_tstep, only.values=T)$values[1])
    
    # country specific timeliness curve
    country_timeliness <- timeliness [country_code == iso3 & !is.na(age), timeliness]
    timeliness_ages    <- timeliness [country_code == iso3 & !is.na(age), age]
    
    # Beta only a single file
    s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
    jt        <- 3 # how many ages to expand to s (or to weeks)
    beta_full <- matrix (0, 
                         ncol = 254, 
                         nrow = 254) # expanding to 100 years; contacts> 80 set to previously 0
    
    beta      <- contact_tstep
    r0_tstep  <- Re (eigen (contact_tstep, only.values=T)$values[1]) # will need it to make sure R0 is the same after rescaling
    
    #create a new contact matrix for each age stratum in model
    
    #expand contact of first 3 age-strata with itself
    
    beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
      A = beta[1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same,
      expand_rows =  s, expand_cols =  s,
      rescale_rows = FALSE, rescale_cols = FALSE
    )
    
    # expand contact for first 3 age-strata with all other contacts
    beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
      A = beta[1:jt,(jt+1):ncol(beta)]/s,
      expand_rows = s, expand_cols = 1, 
      rescale_rows = F, rescale_cols = F)
    
    beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
      A = beta[(jt+1):nrow(beta),1:jt],
      expand_rows = 1, expand_cols = s,
      rescale_rows = F, rescale_cols = F)
    
    beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-  beta[(jt+1):nrow(beta),(jt+1):ncol(beta)]
    
    # make the matrix symmetric
    beta_full <- (beta_full + t(beta_full)) / 2
    #	beta_full_country <- (beta_full + t(beta_full)*(pop.vector_full%*%t(1/pop.vector_full)))/2
    
    # make sure the R0 is what it should be
    beta_full <- (r0_tstep / Re(eigen(beta_full, only.values=T)$values[1])) * beta_full
    
    # Reduce file-writes
    # beta_full <- sweep(beta_full, 2, pop.vector_full, "/") # these are the infection rates from Wallinga et al
    # here, rows should be divided by the appropriate population vector, not columns, but Fortran uses rows, rather than columns to calculate force of infection so it works out.
    
    beta_file <- paste0 ("outcome/",
                         save.scenario, 
                         "/",
                         foldername,
                         "/input/",
                         psadir,
                         fortran_country_code,
                         "_001beta.txt")
    
    # generate DynaMICE input files for each year
    for (y in years) {
      
      #old script groups those aged 70-80, but division is by actual popsize
      pop.vector <- population[country_code == iso3 & year == y, value]
      
      # first expand polymod matrix (contact_tstep) and poopulation vector and
      # then divide by population sizes, otherwise it doesn't work.
      
      pop.vector_full <- c(rep(pop.vector[1:jt]/s, each = s), pop.vector[(jt +1):length(pop.vector)])
      # change zero values to 1 to avoid division by 0 
      pop.vector_full[pop.vector_full==0] <- 1
      
      
      if (vaccination >= 1) {
        
        # Maximum coverage can (obviously) only be 100%
        # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
        # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
        # In essence, this becomes the inverse of the cumulative timeliness curve
        cycov <- coverage_routine [country_code == iso3 & year == y & vaccine == "MCV1", coverage]/timeliness[country_code == iso3 & is.na(age), prop_final_cov]
        
        if (length(cycov) == 0) {   # check if vaccine is not yet introduced and thereby, coverage value missing for this year
          cycov <- 0
        } else if (is.na(cycov)) {
          cycov <- 0
        }
        
        country_year_timeliness_mcv1 <- 1 - min(
          cycov,
          1
        ) * country_timeliness
        
        country_year_timeliness_mcv1 <- (
          country_year_timeliness_mcv1[1:(length(country_year_timeliness_mcv1)-1)] - country_year_timeliness_mcv1[2:(length(country_year_timeliness_mcv1))]
        )/(country_year_timeliness_mcv1[1:(length(country_year_timeliness_mcv1)-1)])
        
        # Timeliness is reported by week in the first year, and by month in the second year. Assume there is no vaccination in between
        country_year_timeliness_mcv1[is.na(country_year_timeliness_mcv1)] <- 0
        country_year_timeliness_mcv1[is.nan(country_year_timeliness_mcv1)] <- 0
        country_year_timeliness_mcv1_allages <- rep(0, 254)
        country_year_timeliness_mcv1_allages[round(timeliness_ages)] <- country_year_timeliness_mcv1
      } else {
        country_year_timeliness_mcv1_allages <- rep(0, 254)
      }
      
      if(vaccination == 2){
        country_year_mcv2 <- coverage_routine [country_code == iso3 & year == y & vaccine == "MCV2", coverage]
      } else {
        country_year_mcv2 <- 0
      }
      
      if ( length (country_year_mcv2) == 0 ) {
        country_year_mcv2 <- 0
      }
      
      # if ( is.na(country_year_mcv2) ){
      #   country_year_mcv2 <- 0
      # }
      
      
      # First three ages are modelled in weekly strata
      pop_year <- c(
        rep (population[country_code == iso3 & year == y & age_to == 0, value]/52, 52),
        rep (population[country_code == iso3 & year == y & age_to == 1, value]/52, 52),
        rep (population[country_code == iso3 & year == y & age_to == 2, value]/52, 52),
        population [country_code == iso3 & year == y & age_to >= 3, value]
      )
      
      # change any zeros to 0 to avoid division by zero
      pop_year [pop_year == 0] <- 1
      
      # write data for each year
      i <- which (years == y)
      
      if (i<10) {
        # beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
        #                    psadir,fortran_country_code,"_00",i,"beta.txt")
        dynamice_input_file <- paste0 ("outcome/",
                                       save.scenario, 
                                       "/",
                                       foldername,
                                       "/input/",
                                       psadir,
                                       fortran_country_code,
                                       "_00",
                                       i,
                                       "measle_data.txt")
      } else if(i<100) {
        # beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
        #                    psadir,fortran_country_code,"_0",i,"beta.txt")
        
        dynamice_input_file<- paste0 ("outcome/",
                                      save.scenario, 
                                      "/",foldername,
                                      "/input/",
                                      psadir,
                                      fortran_country_code,
                                      "_0",
                                      i,
                                      "measle_data.txt")
      } else {
        # beta_file <- paste0("outcome/",save.scenario, "/",foldername,"/input/",
        #                    psadir,fortran_country_code,"_",i,"beta.txt")
        
        dynamice_input_file <- paste0 ("outcome/",
                                       save.scenario, 
                                       "/",
                                       foldername,
                                       "/input/",
                                       psadir,
                                       fortran_country_code,
                                       "_",
                                       i,
                                       "measle_data.txt")
      }
      
      fwrite(
        as.list(
          as.data.table(
            format(
              beta_full,
              digits=14,
              scientific=FALSE
            )
          )
        ),
        beta_file,
        quote=FALSE,
        row.names=FALSE,
        col.names=FALSE,
        sep=" "
      )
      
      dynamice_input_vector <- c(
        t_end,
        1,
        gamma,
        tstep,
        take1,
        take2,
        take3,
        degree1,
        degree2,
        degree3,
        length(sia_a0),
        sia_a0,
        r0_target,
        amplitude,
        1,
        length(sia_a1),
        sia_a1,
        t_end,
        t_cov,
        length(t_sia),
        t_sia,
        length(country_year_timeliness_mcv1_allages),
        country_year_timeliness_mcv1_allages,
        country_year_mcv2,
        length(sia_coverage),
        sia_coverage,
        length(t_run),
        t_run,
        format(
          pop_year,
          digits=14
        )
      )
      
      fwrite (
        list ("a"=dynamice_input_vector),
        file = dynamice_input_file,
        col.names = FALSE
      )
    }
    
    # run actual model
    if (psa == 0) {
      p <- 0
    } else {
      p <- r
    }
    writelog ("gavi_log", paste0(iso3, "; Run ",r,"/",runs,"; Finished generating data"))
    
    # process debugging options
    if (debug_country != "*" & !(iso3 %in% debug_country)){
      debug_debug			    <- 0
      debug_compartments	<- 0
      debug_age		       	<- 0
      debug_timepoints  	<- 0
      debug_relative	  	<- 0
    } else {
      debug_debug 	     	<- as.integer (debug_spinup) + 2 * as.integer(debug_model)
      debug_compartments	<- as.integer (debug_compartments)
      debug_relative		  <- as.integer (debug_relative)
    }
    
    # run the model
    writelog ("gavi_log",paste0(iso3, "; Run ",r,"/",runs,"; Start model"))
    updateProgress (iso3, ii, runs, r, 2)
    
    if(Sys.info()[["sysname"]] == "Windows"){
      model <- paste0(
        'cd "',
        'model/compiled/" & ',
        measles_model, " ",
        length(years), " ",
        tstep, " ",
        save.scenario, " ",
        foldername, " ",
        fortran_country_code, " ", 
        p, " ",
        "WIN", " ",
        debug_debug, " ",
        debug_compartments, " ",
        debug_age, " ",
        debug_timepoints, " ",
        debug_relative
      )
      
      model <- gsub ("/", "\\\\", model)
      shell (model, intern = TRUE)
      
    } else {
      model <- paste0 (
        'cd "', 
        'model/compiled/"; ./',
        measles_model, " ",
        length(years), " ",
        tstep, " ",
        save.scenario, " ",
        foldername, " ",
        fortran_country_code, " ", 
        p, " ",
        "LIN", " ",
        debug_debug, " ",
        debug_compartments, " ",
        debug_age, " ",
        debug_timepoints, " ",
        debug_relative
      )
      system(model, intern=TRUE)
    }
  }
  
  writelog ("gavi_log", paste0(iso3, "; Run ",r,"/",runs,"; Finished model"))
  updateProgress (iso3, ii, runs, r, 3)
  
  # remove input data for this country
  # if files are removed around the same time by different workers, position may change
  if(remove_files){
    files <- list.files(
      paste0("./outcome/",save.scenario, "/",foldername,"/input/",psadir), full.names = TRUE
    )
    do.call(
      file.remove, list(
        files[
          grepl(
            fortran_country_code,
            files
          )
        ]
      )
    )
  }
  
  if (process_results){
    #process model output
    cases <- fread(
      paste0(
        "outcome/",save.scenario, "/",
        foldername,
        "/output/",
        psadir,
        fortran_country_code,
        "_age_stratified_cases_byyear.txt"
      )
    )
    colnames(cases) <- as.character(c(1:ncol(cases))-1)
    cases[,"year"] <- years
    cases <- melt(cases, id.vars = "year", variable.factor=FALSE, value.factor = FALSE, variable.name="age", value.name="cases")
    class(cases$year) <- "numeric"
    class(cases$age) <- "numeric"
    class(cases$cases) <- "numeric"
    
    dat_cfr <- cfr[country_code == iso3, CFR]
    if(psa > 0){
      dat_cfr <- dat_cfr * as.numeric(psa_var [r, c("mortality_input")])
    }
    
    # deaths for those under 5: CFR/100 * cases
    # deaths for those over 5 and under 10: CFR/200 * cases
    # deaths for those over 10: 0
    cases[, "deaths"] <- 0
    cases[age < 10, "deaths"] <- cases[age < 10, cases] * (dat_cfr/200)
    cases[age < 5, "deaths"] <- cases[age < 5, cases] * (dat_cfr/100)
    
    for(y in years){
      dat_life_exp <- lexp[country_code == iso3 & year == y, value] # DALYs CHANGE
      z <- 1
      while(length(dat_life_exp) == 0){
        dat_life_exp <- lexp[country_code == iso3 & year == years[which(years==y) - z], value]
        z <- z+1
        if(z == length(years)){
          break
        }
      }
      cases[year == y, "cohort_size"] <- population[country_code == iso3 & year == y, value]
      
      # Change 1 - Use life expectancy that changes by year (old lexp was fixed at life expectanct of 2017)
      
      # Change 2 - 
      # In the old version, assume that average age of the deaths was 4 years, and use lexp - 4 to estimate years of life lost for each death
      # In the new version, use the actual age of the child that died and life-expectancy of that year to estimate the effect
      
      # cases[year == y, "dalys"] <- (cases[year == y, cases] - cases[year == y, deaths]) * 0.002 + cases[year == y, deaths] * (dat_life_exp - cases[year == y, age])
      cases[year == y, "dalys"] <- (cases[year == y, cases] - cases[year == y, deaths]) * 0.002 + cases[year == y, deaths] * (dat_life_exp - 4) # DALYs CHANGE
    }
    
    cases[, "country_code"] <- iso3
    cases <- cases[, c("year", "age", "country_code", "cohort_size", "deaths", "cases", "dalys"), with=F]
    
    setorder (cases, year, age)
    
    writelog ("gavi_log", paste0 (iso3, "; Run ", r, "/", runs, "; Return results"))
    updateProgress (iso3, ii, runs, r, 4)
    
    return (cases)
  }
  
} # end of function -- runCountry
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# runScenario
#
# Run the measles model for a given scenario
# ------------------------------------------------------------------------------
runScenario <- function (vaccine_coverage_folder    = "", 
                         coverage_prefix            = "",
                         touchstone                 = "",
                         antigen                    = "",
                         scenario_name,
                         scenario_number,
                         vaccine_coverage_subfolder = "",
                         burden_template,                   # burden template file
                         burden_estimate_folder,            # burden estimate folder  
                         group_name,                        # modelling group name
                         countries                  = "all", 
                         cluster_cores              = 1,
                         psa                        = 0,    # psa runs; 0 for single run
                         vaccination,  # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
                         using_sia,    # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA
                         measles_model,                     # measles model
                         debug_model                = FALSE # debug model (T/F)
) {
  
  # changes 2019: 
  #       a) add scenario folder pathname for clarity and for easier running on the cluster (fortran code changed too)
  #       b) change contact matrix to fix augumented mixing in <3 yo and to reflect country specific mixing (rescale polymod by local population structure to keep contacts reciprocal)
  #          - uses newly rescaled contact matrix to expand polymod to yearly ages up to 100
  #       c) time-varying CFRs (postprocessing)
  #       d) MCV1 and MCV2 dependency - give MCV2 to those who received MCV1 for MCV2 coverage <= MCV1, if MCV2 coverage > MCV1 distribute the remainder randomly
  #       e) MCV1 and SIA dependency - use linear regression of data from Portnoy et al 2018 to inform the proportion of zero-dose children reached by SIA campaign
  #       f) due to issues with age_range_verbatim field in the scenario (coverage) input files, these runs use age_first and age_last
  
  # save to the scenario folder   
  save.scenario <- scenario_number
  
  ages 		    <- c (0:100)	              # Numeric vector with age-strata that are reported (Model ALWAYS models 101 age-groups; 0-2yo in weekly age-strata, and 3-100yo in annual age-strata)
  psa 		    <- psa		                  # Number of PSAs to run (0 to use deterministic).
  
  # --------------------------------------------------------------------------
  # input data files
  #
  # vaccine coverage file for routine immunisation
  data_coverage_routine <- paste0 (vaccine_coverage_folder, 
                                   vaccine_coverage_subfolder,
                                   "routine_",
                                   scenario_name, 
                                   ".csv")
  
  # separate vaccine coverage file for SIA (supplementary immunisation activities)
  data_coverage_sia <- paste0 (vaccine_coverage_folder, 
                               vaccine_coverage_subfolder,
                               "sia_",
                               scenario_name, 
                               ".csv")
  
  # file with timeliness data (csv; in ./input/)
  data_timeliness <- "timeliness_new.csv"
  
  # file with R0 (csv; in ./input/)
  data_r0 		<- "r0.csv" #"r0est.csv"
  
  # file with life-expectancies (csv; in ./input/)
  data_life_exp 	<- "demographicdata2019/201910gavi-4_dds-201910_lx0_both.csv"
  # data_life_exp 	<- "new/lx_old.csv"
  
  # template Gavi (used to check if all data is present)
  # data_gavi_template <- "central-burden-template_2019104.csv" 
  data_gavi_template <- burden_template
  
  # file with CFRs (csv; in ./input/)
  data_cfr 		<- "cfrs_new_noage.csv"
  
  # file with population sizes (csv; in ./input/)
  data_pop 		<- "demographicdata2019/201910gavi-4_dds-201910_2_int_pop_both.csv"
  #data_pop 		<- "new/pop_old.csv"
  
  # file with contact_data
  data_contact	<- "contact/uk_polymod_physical_101.csv"  # use newly rescaled polymod to yearly age bands to 100
  
  # data with PSA variables
  # should be the same for each scenario. Will be CREATED if does not exist
  data_psa 		<- "psa_variables.csv"
  
  # expected remaining years of life
  # data_life_exp_remain <- "demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both.csv"
  data_life_exp_remain <- "demographicdata2019/201910gavi-4_dds-201910_2_life_ex_both_full.csv"
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # Advanced options 
  # --------------------------------------------------------------------------
  
  sia.method	<- 1				             # 1 = variable take; 2 = variable degree
  dinf		    <- 14				             # duration of infection (days)
  amplitude 	<- 0.05			             # amplitude for seasonality
  
  # take 		    <- c (0.85, 0.95, 0.98)  # vaccine efficacy for take1 (take dose 1, before age 1), take2 (dose 1, after age 1) & take3 (dose 2). Note that dose2 only has an effect if vaccine==2.
  # used to mean before: take [1] # vaccine efficacy (dose 1, before age 1)
  #                      take [2] # vaccine efficacy (dose 1, after age 1)
  # updated meaning:     take [1] refers to vaceffbyage_a (intercept)
  #                      take [2] refers to vaceffbyage_b (slope)
  take 		    <- c (0.64598, 0.01485, 0.98)
  
  degree 		  <- c (0.85, 0.95, 0.98)  # vaccine efficacy for degree1 (degree dose 1, before age 1), degree2 (dose 1, after age 1) & degree3 (dose 2). Note that dose2 only has an effect if vaccine==2.			
  tstep			  <- 1000				           # Number of time steps in a year
  
  
  # ----------------------------------------------------------------------------
  # Measles model
  # filename of compiled fortran-model (in ./model/compiled/)
  # ----------------------------------------------------------------------------
  # measles_model <- "vaccine2019_sia_singlematrix" # change this to reflect the right version of the fortran code
  # measles_model <- "vaccine2019_sia_singlematrix.exe" # change this to reflect the right version of the fortran code
  measles_model <- measles_model
  # ----------------------------------------------------------------------------
  
  # number of clusters to use
  # if larger than 1, country-specific model runs are distributed over specified number of clusters
  # note that model uses a lot of memory, so might not want to max out all clusters
  use_cluster  <- cluster_cores   # debug #
  remove_files <- TRUE
  
  # may want to process results after generating all data. Note OUTPUT files are not removed if remove_files == TRUE and process_results == FALSE
  process_results <- F
  run_model       <- TRUE
  # folder will be created if not given - should usually be commented out, except when run model is FALSE
  
  
  # ------------------------------------------------------------------------------
  # Debug
  # ------------------------------------------------------------------------------
  
  debug_country		  <- "*"			#ISO3 codes of country to debug, * to debug all countries
  debug_spinup		  <- FALSE		#TRUE/FALSE: If true, generate data for spin-up period of model
  
  debug_model       <- debug_model
  # debug_model			  <- FALSE		#TRUE/FALSE: If true: generate data for period after spin-up
  
  debug_compartments<- 1
  # debug_compartments<- 0			  #TRUE/FALSE: If true: output size of each compartment. If false: output number of cases. If 2: debug vaccinated
  debug_age         <- 0        #0-2. If 0: output all in annual age-strata. If 1: output age 0-2 in weekly age-strata, 3-100 in annual age-strata. If 2: sum all age-strata.
  debug_timepoints	<- 0			  #0-2. If 0: output per year. If 1: output per timepoint and report first 25% of timepoints. If 2: output per timepoint and report all timepoints.
  debug_relative		<- FALSE		#If true: output proportion of new cases. If false, output absolute number of new cases.
  
  # START OF MODEL  
  # ----------------------------------------------------------------------------
  # SETUP
  
  # load correct libraries when on cluster (using open MPI)
  if ("cluster_using_openmpi" %in% commandArgs()){
    using_openmpi <- TRUE
    require("doMPI")
    require("parallel")
    cl <- startMPIcluster()
    registerDoMPI(cl)
    print("using openMPI")
  } else {
    using_openmpi <- FALSE
    print("no openMPI")
  }
  
  # clean environment
  if (file.exists ("gavi_progress")) {
    file.remove ("gavi_progress")
  }
  
  # are we running a PSA - a deterministic or "stochastic" model?		
  if (psa == 0) {
    det_stoch <- "deter" # deterministic
    runs      <- 1
    
  } else {
    det_stoch <- "stoch" # stochastic
    runs      <- psa
  }
  
  # create folder with correct name for in- and output data if not yet exists
  # typically foldername should not exist - but may come in handy when only processing results
  if ( !exists("foldername_analysis") ) {
    
    foldername <- paste0 (
      format(Sys.time(),format="%Y%m%d"),
      "_v",
      vaccination,
      "_s",
      using_sia,
      "_",
      det_stoch
    )
    
    dir.create(
      file.path(
        paste0(
          getwd(),
          "/outcome/", save.scenario, "/",
          foldername
        )
      ), recursive = T
    )
    
    dir.create(
      file.path(
        paste0(
          getwd(),
          "/outcome/", save.scenario, "/",
          foldername,
          "/input"
        )
      ), recursive = T
    )
    
    dir.create(
      file.path(
        paste0(
          getwd(),
          "/outcome/", save.scenario, "/",
          foldername,
          "/output"
        )
      ), recursive = T
    )
    
    if (psa > 0) {
      
      for(p in 1:psa){
        if(p<10){
          p <- paste0("00",p)
        } else if(p<100){
          p <- paste0("0",p)
        }
        
        dir.create(
          file.path(
            paste0(
              getwd(),
              "/outcome/", save.scenario, "/",
              foldername,
              "/input/",
              paste0("run",p)
            )
          ), recursive = T
        )
        
        dir.create(
          file.path(
            paste0(
              getwd(),
              "/outcome/", save.scenario, "/",
              foldername,
              "/output/",
              paste0("run",p)
            )
          ), recursive = T
        )
      }
    }
  } else {
    foldername <- foldername_analysis
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  
  # log
  writelog ("gavi_log",paste0("Main; gavi.r started"))
  
  if (using_openmpi) {
    writelog ("gavi_log", paste0 ("Main; OpenMPI enabled"))
  } else {
    writelog ("gavi_log", paste0 ("Main; OpenMPI disabled"))
  }
  
  # read_data
  coverage_routine	<- fread (paste0 (data_coverage_routine))
  coverage_sia		  <- fread (paste0 (data_coverage_sia))
  timeliness			  <- fread (paste0 ("input/", data_timeliness))
  cfr					      <- fread (paste0 ("input/", data_cfr))
  rnought			    	<- fread (paste0 ("input/", data_r0))
  population		  	<- fread (paste0 ("input/", data_pop))
  lexp 		      		<- fread (paste0 ("input/", data_life_exp))
  template		    	<- fread (paste0 (data_gavi_template))
  contact			    	<- fread (paste0 ("input/", data_contact))
  lexp_remain       <- fread (paste0 ("input/", data_life_exp_remain)) 
  
  # coverage_sia has multiple entries per year for some countries, take only the first entry
  # coverage_sia <- coverage_sia[, .SD[1], by = c("country_code", "year")]
  
  # --------------------------------------------------------------------------
  # if psa variables file does not exist, then create psa variables file
  if (psa > 0) {
    if (file.exists (paste0 ("input/", data_psa))) {
      
      # read csv if file already exists
      psa_var <- fread (paste0 ("input/", data_psa))
      
      # check if psa_var corresponds with psa
      if(nrow(psa_var) != psa){
        stop(paste0("Number of runs in PSA file is not the same as those specified! Variables used in the same run in each scenario should be similar. Delete or rename './input/",data_psa," if a new file needs to be created, or change the number of PSAs to ",nrow(psa_var)," if the same file should be used."))
      }
    } else {
      
      # create csv if file does not exist
      psa_var <- data.table(
        run_id = c(1:psa),
        take1_input     = rep (NA, psa),
        take2_input     = rep (NA, psa),
        take3_input     = rep (NA, psa),
        mortality_input = rep (NA, psa)
      )
      
      for(t in 1:3) {
        
        # use 5% difference in take
        tk <- runif (n   = psa, 
                     min = (take[t] * 0.95), 
                     max = (take[t] * 1.05))
        
        # vaccine efficacy is bounded between 0 and 1
        tk [which(tk < 0)] <- 0
        tk [which(tk > 1)] <- 1
        
        psa_var [, paste0 ("take", t, "_input")] <- tk
      }
      
      # use 25% difference in mortality, will be multiplied by CFR in each country
      psa_var [, "mortality_input"] <- runif (n  = psa,
                                             min = 0.75,
                                             max = 1.25)
      
      fwrite (psa_var, paste0 ("input/", data_psa))
    }
  }
  # --------------------------------------------------------------------------
  
  # process data
  
  # if countries are specified to all, then set countries to all countries in coverage file
  if (countries == "all") {
    countries	<- as.character (unique (coverage_routine [, country_code] ) )  
  }
  
  # start and end years (should go in as input to function -- INPUT-FUNCTION)
  years <- as.numeric (c(1980:2100))
  
  for(c in countries) {
    if(psa > 1){
      for(p in 1:psa){
        write(
          paste0(
            c,
            " ",
            paste0(c(rep(0,(nchar(psa)-nchar(p))),p),collapse=""),
            " 0"
          ),
          file   = "gavi_progress",
          append = TRUE
        )
      }
    } else {
      write(
        paste0(
          c,
          " ",
          "1 0"
        ),
        file   = "gavi_progress",
        append = TRUE
      )
    }
  }
  
  contact <- contact [ , -"contact.age"]
  contact <- as.matrix (contact)
  
  
  r0_basic <- Re (eigen (contact * dinf, only.values = T)$values[1])
  gamma    <- 1 / (dinf * tstep/365)
  
  # Run model
  writelog ("gavi_log", paste0 ("Main; Foldername: ", foldername))
  
  if(using_openmpi | use_cluster > 1){
    if(!using_openmpi){
      
      if ( use_cluster != detectCores() ) {
        warning (paste0(detectCores(), " cores detected but ", use_cluster, " specified."))
      }
      
      cl <- makeCluster (use_cluster)
      registerDoParallel (cl)
      
    } else {
      print(paste0("Clustersize: ", clusterSize(cl)))
    }
  }
  
  # foreach will run countries and PSA runs in parralel if a parallel backend 
  # is registered, and sequentially otherwise
  combine <- foreach (
    ii = 1:length(countries),
    .packages = c("data.table"),
    .errorhandling="pass", 
    .export = c("runCountry", "writelog", "expandMatrix", "updateProgress")
  ) %:% foreach (
    r = 1:runs,
    .errorhandling = "pass"
  ) %dopar% {
    
    out_run <- runCountry (ii                 = ii, 
                           countries          = countries, 
                           years              = years,
                           vaccination        = vaccination, 
                           using_sia          = using_sia, 
                           dinf               = dinf, 
                           gamma              = gamma, 
                           r0_basic           = r0_basic, 
                           amplitude          = amplitude, 
                           take               = take, 
                           degree             = degree, 
                           sia.method         = sia.method,
                           coverage_routine   = coverage_routine, 
                           coverage_sia       = coverage_sia, 
                           timeliness         = timeliness, 
                           contact            = contact, 
                           rnought            = rnought, 
                           population         = population, 
                           lexp               = lexp, 
                           cfr                = cfr,
                           tstep              = tstep,
                           save.scenario      = save.scenario, 
                           foldername         = foldername,  
                           measles_model      = measles_model, 
                           debug_model        = debug_model, 
                           debug_spinup       = debug_spinup, 
                           debug_age          = debug_age, 
                           debug_compartments = debug_compartments, 
                           debug_country      = debug_country, 
                           debug_relative     = debug_relative, 
                           debug_timepoints   = debug_timepoints,
                           r                  = r,
                           runs               = runs, 
                           psa                = psa, 
                           psa_var            = psa_var,
                           process_results    = process_results, 
                           run_model          = run_model, 
                           remove_files       = remove_files
    )
    return(out_run)
  }
  
  if(using_openmpi | use_cluster > 1){
    if(!using_openmpi){
      stopCluster(cl)
    } else {
      closeCluster(cl)
    }
  }
  
  # check for errors
  errorcount <- 0
  
  for (i in 1:length(combine)){
    if("error" %in% class(combine[[i]])){
      errormessage <- paste0("Error in task ",i,": ",combine[[i]])
      warning(errormessage)
      writelog(paste0(gavi.dir,"gavi_log"),errormessage)
      errorcount <- errorcount + 1
      #remove from data
      combine[[i]] <- NULL
    }
  }
  #if(errorcount > 0){
  #	stop(paste0("There were ",errorcount," errors."))
  #}
  
  for(ii in 1:length(countries)){
    for(r in 1:runs){
      if(ii==1 & r==1){
        out_gavi_agestratified <- combine[[1]][[1]]
      } else {
        out_gavi_agestratified <- rbindlist(
          list(
            out_gavi_agestratified,
            combine[[ii]][[r]]
          )
        )
      }
    }
  }
  
  if (process_results){
    #clean environment
    if(remove_files){
      if(psa >0){
        for(p in 1:psa){
          if(p<10){
            p <- paste0("00",p)
          } else if(p<100){
            p <- paste0("0",p)
          }
          do.call(file.remove, list(list.files(paste0("./outcome/", save.scenario, "/",foldername,"/output/run",p,"/"), full.names = TRUE)))
        }
      } else {
        do.call(file.remove, list(list.files(paste0("./outcome/", save.scenario, "/",foldername,"/output/"), full.names = TRUE)))
      }
      
      if(Sys.info()[["sysname"]] == "Windows"){
        do.call(file.remove, list("./model/compiled/fort.6"))
      }
    }
    
    # write output
    if(psa >0){
      det_stoch <- "stochastic"
    } else {
      det_stoch <- "deterministic"
    }
    
    # make output similar to template
    # only report years to be reported
    out_gavi_agestratified <- out_gavi_agestratified[year %in% unique(template[, year])]
    out_gavi_agestratified[, country_name := unique(template[country == country_code, country_name]), by = country_code]
    out_gavi_agestratified[, "disease"] <- "Measles"
    colnames(out_gavi_agestratified)[3] <- "country"
    out_gavi_agestratified <- out_gavi_agestratified[, colnames(template), with=F]
    setorder(out_gavi_agestratified, country, year, age)
    
    # write data
    filename <- paste0(
      "outcome/", save.scenario, "/",
      foldername,
      "/",
      format(Sys.time(),format="%Y%m%d"),
      "_gavi_measles_",
      c("vax_none","vax_mcv1","vax_mcv1_mcv2")[vaccination+1],
      c("_sia_no","_sia_yes")[using_sia+1],
      "_",
      det_stoch,
      "_agestratified.csv"
    )
    
    z <- 1
    
    while (file.exists(filename)) {
      warning(paste0(filename)," already exists, adding _",z)
      
      filename <- paste0(
        "outcome/", save.scenario, "/",
        foldername,
        "/",
        format(Sys.time(),format="%Y%m%d"),
        "_gavi_measles_",
        c("vax_none","vax_mcv1","vax_mcv1_mcv2")[vaccination+1],
        c("_sia_no","_sia_yes")[using_sia+1],
        "_",
        det_stoch,
        "_agestratified_",
        z,
        ".csv"
      )
      
      if(!file.exists(filename)){
        break
      }
      z <- z + 1
    }
    
    fwrite(
      out_gavi_agestratified,
      filename,
      row.names=FALSE
    )
  }
  
  # ----------------------------------------------------------------------------
  # psa data (and psa variables file) are generated once at the start of the program
  #
  # the following code is not needed -- commented section
  # write file stochastic variables
  # if(psa > 0){
  #   stoch_file <- paste0(format(Sys.time(),format="%Y%m"),"_out_measles_stochastic_variables.csv")
  #   if(file.exists(paste0("outcome/",stoch_file))){
  #     #read csv if file already exists
  #     stoch_file <- read.csv(paste0("input/",data_psa))
  #     #check if psa_var corresponds with psa
  #     if(nrow(stoch_file) != psa){
  #       stop(paste0("Number of runs in PSA output file is not the same as those specified! File is NOT overwritten, please delete or rename the old file if a new file needs to be generated"))
  #       writelog("gavi_log",paste0("Main; gavi.r aborted PSA error"))
  #       if(using_openmpi){
  #         mpi.quit()
  #       }
  #     }
  #   } else {
  #     # create csv if file does not exist
  #     print("Creating stochastic parameters file")
  #     stoch_file <- psa_var[,c("run_id","take1_input","take2_input","take3_input")]
  #     #get country specific CFR
  #     mortality <- as.numeric(psa_var[,"mortality_input"])
  #     for(c in 1:length(countries)){
  #       cfr <- Crit[which(Crit$country==as.character(countries[c])),"CFR"]
  #       stoch_file[,paste0(countries[c],"_CFR")] <- cfr*mortality
  #     }
  #     fwrite(
  #       stoch_file,
  #       paste0(
  #         "outcome/",
  #         format(Sys.time(),format="%Y%m"),
  #         "_out_measles_stochastic_variables.csv"
  #       ),
  #       row.names=FALSE
  #     )
  #   }
  # }
  # ----------------------------------------------------------------------------
  
  
  # process results
  
  report_years <- sort(unique(template$year))
  ages         <- sort(unique(template$age))
  
  # get years that the model was run for (can be different from reporting years)
  years 		   <- sort(unique(as.numeric(coverage_routine[,year])))
  
  
  # set cfrs before 2000 to 2000
  cfr.year <- rbindlist(lapply(1980:1999, function(i) copy(cfr[Year ==2000,])[,Year := i]))
  cfr <- rbind(cfr.year, cfr)
  
  # set Kosovo mortaity to Serbia as it was once Serbia
  cfr.xk <- rbindlist(lapply("XK", function(i) copy(cfr[Code == "SRB",])[,Code := i]))
  cfr.xk[, Country := "Kosovo"]
  cfr <- rbind(cfr, cfr.xk)
  
  # expand cfr by age so that it uses the right value for <5 and 5-9 and 0 for >=10
  cfr.year.all <- rbindlist(lapply(0:100, function(i) copy(cfr)[, age:=i]))
  over10 <- T
  
  if (over10){
    cfr.year.all[, cfr.value :=  if (age < 5) under5 else if (age > 4 & age <10) over5 else 0, 
                 by = c("Code", "Year", "age") ] 
  } else {
    cfr.year.all[, cfr.value :=  if (age < 5) under5 else if (age > 4 ) over5, 
                 by = c("Code", "Year", "age") ]
  }
  cfr.year.all$Year = as.integer(cfr.year.all$Year)
  
  
  # ANALYSE RUNS 
  
  # get the names of all the files in all the subfolders the stochastic runs specified in scenarioname; separate cases and popsize files
  # will use those to read them all in at once and put in a data.table
  myfiles <- list.files(path = paste0("outcome/", save.scenario, "/",foldername,"/output/"), 
                        recursive = T, pattern = "cases", full.names = T)
  myfiles.popsize <- list.files(path = paste0("outcome/", save.scenario, "/",foldername,"/output"), 
                                recursive = T, pattern = "popsize", full.names = T)
  
  # read in all those files specified in myfiles and put them in a single data table
  all_cases  <- rbindlist(lapply(myfiles, function(fn, ...) {
    res <- withCallingHandlers(
      fread(fn, stringsAsFactors=F, check.names = F, fill = T, 
            col.names = as.character(c(0:100))),
      warning = function(w) { warning(w, fn); }
    )
    res[, country := gsub("^.+/(\\w+)_age.+$","\\1", fn) ]             # get the country code from the filename
    # get run_id for psa runs
    if (psa > 0) {
      res[, run_id := as.integer(gsub("^.+run(\\d{3}).+$","\\1", fn)) ]  # get the run_id from the filename (subfolder)
    }
    res[, year := years]                                               # add year of simulation (from coverage data file)
  })) 
  
  # switch back Kosovo to XK otherwise montagu returns error on upload
  all_cases[country == "XKX", country := "XK"]
  
  all_popsize  <- rbindlist(lapply(myfiles.popsize, function(fn, ...) {
    res <- fread(fn, stringsAsFactors=F, check.names = F, col.names = as.character(c(0:100)))
    res[, country := gsub("^.+/(\\w+)_age.+$","\\1", fn) ]             # get the country code from the filename
    if (psa > 0) {
      res[, run_id := as.integer(gsub("^.+run(\\d{3}).+$","\\1", fn)) ]  # get the run_id from the filename (subfolder)
    }
    res[, year := years]                                               # add year of model was run for (from coverage data file)
  })) 
  all_popsize[country == "XKX", country := "XK"]
  
  if (psa > 0) {
    column_names <- c("year", "country", "run_id")
  } else {
    column_names <- c("year", "country")
  }
  
  # all_cases.m <- melt(all_cases, id.vars = c("year","country"), 
  all_cases.m <- melt(all_cases, id.vars = column_names,
                      measure.vars = c(as.character(0:100)), 
                      variable.name = "age", value.name = "cases", variable.factor = F)
  
  # all_popsize.m <- melt(all_popsize, id.vars = c("year", "country"),
  all_popsize.m <- melt(all_popsize, id.vars = column_names,
                        measure.vars = c(as.character(0:100)), 
                        variable.name = "age", value.name = "cohort_size", variable.factor = F)
  
  # melt produces character variable when variable.factor is set to FALSE - change it to integer
  all_cases.m   [, age := lapply(.SD, as.integer), .SDcols = "age"]
  all_popsize.m [, age := lapply(.SD, as.integer), .SDcols = "age"]
  
  # merge cases and cohort sizes
  # all_runs = merge (all_cases.m, all_popsize.m, by = c("year", "age", "country"), all.x = T)
  all_runs = merge (all_cases.m, all_popsize.m, by = c(column_names, "age"), all.x = T)
  
  # add country_name, life expectancy and disease (Measles) to match template file
  country_names <- unique (subset (template, select = c("country", "country_name")))
  c_names <- country_names$country_name; names(c_names) = country_names$country
  life.exp2 <- rbindlist(lapply(2100, function(i) copy(lexp[year == 2099])[, year:=i]))
  life.exp <- rbind(lexp, life.exp2)
  life.exp.all <- rbindlist(lapply(0:100, function(i) copy(life.exp)[, age:=i]))
  colnames(life.exp)[8] <- "LE"
  # all_runs[, c("country_name", "disease", "LE") := list(c_names[country], "Measles", life.exp[country])]
  
  all_runs[, c("country_name", "disease") := list(c_names[country],"Measles")]
  
  
  # merge mortality rates from CFR file
  all_runs <- merge(all_runs, subset(cfr.year.all, 
                                     select = c( "Code", "Year", "age", "cfr.value")),
                    by.x = c("country", "year", "age" ), by.y = c("Code", "Year", "age"), 
                    sort= F, all.x = T)
  
  all_runs <- merge (all_runs, 
                     subset (life.exp, select= c("country_code", "year", "LE")), 
                     by.x = c("country", "year"), 
                     by.y = c("country_code", "year"))
  
  # ----------------------------------------------------------------------------
  # add data column for remaining life expectancy
  if (psa > 0) {
    all_runs <- lexp_remain [all_runs,
                             .(run_id, i.country, year, age, cases, cohort_size, country_name, disease, cfr.value, LE, value),
                             on = .(country_code = country,
                                    age_from    <= age,
                                    age_to      >= age,
                                    year         = year) ]
  } else {
    all_runs <- lexp_remain [all_runs,
                             .(i.country, year, age, cases, cohort_size, country_name, disease, cfr.value, LE, value),
                             on = .(country_code = country,
                                    age_from    <= age,
                                    age_to      >= age,
                                    year         = year) ]
  }
  
  # rename column "i.country" to "country"
  setnames (x = all_runs, old = "i.country", new = "country")
  
  # save a copy of remaining life expectancy values in disease column
  # rename country column
  setnames (all_runs, 
            old = c ("value"),  
            new = c ("remain_lexp") 
  )
  
  # all_runs [, disease := value]
  # all_runs [, value   := NULL]
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # MCV1 coverage
  coverage_routine_MCV1 <- coverage_routine [vaccine == "MCV1"]
  
  # add MCV1 column
  if (psa > 0) {
    all_runs <- coverage_routine_MCV1 [all_runs, 
                                       .(run_id, i.country, i.year, age, cases, cohort_size, country_name, disease, cfr.value, LE, coverage, remain_lexp),
                                       on = .(country_code = country,
                                              year         = year) ]
  } else {
    all_runs <- coverage_routine_MCV1 [all_runs, 
                                       .(i.country, i.year, age, cases, cohort_size, country_name, disease, cfr.value, LE, coverage, remain_lexp),
                                       on = .(country_code = country,
                                              year         = year) ]
  }
  
  # rename column "coverage" to "MCV1"
  setnames (x = all_runs, 
            old = c("i.country", "i.year", "coverage"), 
            new = c("country",   "year",   "MCV1"))
  # ----------------------------------------------------------------------------
  
  # calculate deaths and dalys
  all_runs[, deaths := cases * cfr.value]
  #all_runs[, dalys := (cases - deaths)*0.002 + deaths*(LE - 4)] # assumes all deaths are age 4, 
  all_runs [, dalys := (cases - deaths) * 0.002 + deaths * (LE - age)] # to account for actual age of death
  
  # ------------------------------------------------------------------------------
  # YLL calculation above is fine since most measles burden < 5 years
  # remaining life expectancy estimation can be improved by using remaining life
  # years by age, year and country
  # ------------------------------------------------------------------------------
  
  # OUTPUT RUNS
  
  # don't need all of these columns for VIMC, save only ones that are needed
  if (psa > 0) {
    save.cols <- c("run_id", colnames(template))
  } else {
    save.cols <- c(colnames(template))
  }
  
  # ----------------------------------------------------------------------------
  save.cols <- c(save.cols, "MCV1", "remain_lexp")
  # ----------------------------------------------------------------------------
  
  # keep all columns
  output_runs <- subset(all_runs, year %in% report_years, select = save.cols)
  
  # burden estimate type -- central or stochastic
  if (psa > 0) {
    burden_estimate_type <- "stochastic_burden_estimate_"
  } else {
    burden_estimate_type <- "central_burden_estimate_"
  }
  
  # burden estimate filename
  burden_estimate_file <- paste0 (burden_estimate_type,
                                  antigen, 
                                  group_name, 
                                  scenario_name, 
                                  ".csv")
  
  # save burden estimates to file
  fwrite (x    = output_runs [order(country, year, age)], 
          file = paste0 (burden_estimate_folder, 
                         burden_estimate_file) )
  
  # clean environment
  writelog ("gavi_log", paste0 ("Main; gavi.r finished"))
  
  if (using_openmpi) {
    mpi.quit()
  }
  
  # return burden estimate filename (cases)
  return (burden_estimate_file)
  
} # end of function -- runScenario
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: estimateDeathsDalys 
#
# Estimate deaths by applying CFR (case fatality rates) to estimated cases and
# also calculate DALYs
# ------------------------------------------------------------------------------
# 
#   save results in corresponding cfr_option subfolder
#   append cfr_option to burden estimates file
# ------------------------------------------------------------------------------
estimateDeathsDalys <- function (cfr_option,
                                 burden_estimate_file, 
                                 burden_estimate_folder, 
                                 vimc_scenario, 
                                 portnoy_scenario, 
                                 psa = 0) {
  
  # read burden estimates (primarily cases)
  burden <- fread (file = paste0 (burden_estimate_folder, 
                                  burden_estimate_file) )
  
  # ----------------------------------------------------------------------------
  # use CFRs -- Wolfson
  if (cfr_option == "Wolfson") {
    
    # read CFRs (in percentage) -- Wolfson
    # cfr estimates for ages below 5 years
    # cfrs for ages above 5 years are half the cfr values for ages below 5 years
    cfr	<- fread ("input/cfr.csv")
    
    # change CFR percentages to rates between 0 and 1
    cfr [, CFR := CFR/100]
    
    # add CFR data column to burden estimates
    if (psa > 0) {
      burden <- cfr [burden, 
                     .(run_id, disease, year, age, i.country, country_name, cohort_size, cases, dalys, deaths, CFR, remain_lexp), 
                     on = .(country_code = country) ]
    } else {
      burden <- cfr [burden, 
                     .(disease, year, age, i.country, country_name, cohort_size, cases, dalys, deaths, CFR, remain_lexp), 
                     on = .(country_code = country) ]
    }
    
    # estimate deaths
    # cfrs for ages above 5 years are half the cfr values for ages below 5 years
    burden [age <  5, deaths := cases * CFR  ]
    burden [age >= 5, deaths := cases * CFR/2]
    
    # rename country column
    setnames (burden, 
              old = c ("i.country"),  
              new = c ("country") 
    )
    
    # drop CFR column
    burden [, CFR  := NULL]
  }
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # use CFRs -- Portnoy
  
  if (cfr_option == "Portnoy") {
    
    # read CFRs (rates between 0 and 1) -- Portnoy
    cfr	<- fread ("input/cfr_scenarios.csv")
    # cfr	<- fread ("input/cfrs_new_noage.csv")
    
    # --------------------------------------------------------------------------
    # cfr estimates are for 2000 to 2030
    # if cfr estimates are required for years below or above this range, then 
    # for years below 2000, set cfr estimates of year 2000
    # for years above 2030, set cfr estimates of year 2030
    
    # find minimum and maximum year
    min_year = min (burden [, year])
    max_year = max (burden [, year])
    
    # set cfrs for years before 2000 to cfr values of year 2000
    if (min_year < 2000) {
      cfr.year <- rbindlist (lapply (min_year:1999, function(i) copy (cfr [year == 2000, ])[, year := i]))
      cfr      <- rbind     (cfr.year, cfr, use.names = TRUE)
    }
    
    # set cfrs after 2030 to 2030
    if (max_year > 2030) {
      cfr.year <- rbindlist (lapply (2031:max_year, function(i) copy (cfr [year == 2030, ])[, year := i]))
      cfr      <- rbind     (cfr, cfr.year, use.names = TRUE)
    }
    # --------------------------------------------------------------------------
    
    # rename columns -- cfr of vimc_scenario and portnoy_scenario
    # cfrs for under 5 (< 5) and over 5 (>= 5) years
    setnames (x = cfr, 
              old = c (paste0 ("cfr_under5_", vimc_scenario, "_", portnoy_scenario),      
                       paste0 ("cfr_over5_" , vimc_scenario, "_", portnoy_scenario) ), 
              new = c ("under5", 
                       "over5" ) 
    )
    
    # add CFR data column to burden estimates
    if (psa > 0) {
      burden <- cfr [burden, 
                     .(run_id, disease, year, age, country, country_name, cohort_size, cases, dalys, deaths, over5, under5, remain_lexp), 
                     on = .(country_code = country, 
                            year         = year) ]
    } else {
      burden <- cfr [burden, 
                     .(disease, year, age, country, country_name, cohort_size, cases, dalys, deaths, over5, under5, remain_lexp), 
                     on = .(country_code = country, 
                            year         = year) ]
    }

    # estimate deaths for ages under 5 years
    burden [age < 5, deaths := cases * under5]
    
    # estimate deaths for ages over 5 years
    burden [age >= 5, deaths := cases * over5]
    
    # drop cfr columns -- under5 and over5
    burden [, ':=' (under5 = NULL, 
                    over5  = NULL)]
  }
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # DALYs
  # calculate dalys = (ylds) + (ylls)
  burden [, dalys := ((cases - deaths) * 0.002) + (deaths * remain_lexp)]
  
  # drop column -- remaining life expectancy
  burden [, remain_lexp  := NULL]
  # ----------------------------------------------------------------------------
  
  # append/suffix cfr_option to the end of filename
  updated_burden_estimate_file <- str_replace (string      = burden_estimate_file, 
                                               pattern     = ".csv", 
                                               replacement = paste0 ("_", cfr_option, ".csv")  
  )
  
  # save updated burden estimate file (cases + deaths) to file
  # cfr_option is also the name of the subfolder
  fwrite (x    = burden, 
          file = paste0 (burden_estimate_folder,
                         cfr_option, "/", 
                         updated_burden_estimate_file)
  )
  
  return ()
  
} # end of function -- estimateDeathsDalys
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# diagnostic_plots
#
# diagnostic plots of vaccine coverage and burden estimates (cases, deaths, dalys)
# ------------------------------------------------------------------------------
diagnostic_plots <- function (vaccine_coverage_folder,
                              coverage_prefix,
                              touchstone,
                              antigen,
                              scenarios,
                              base_scenario,
                              burden_estimate_folder,
                              plot_folder,
                              group_name,
                              countries,
                              cfr_options, 
                              psa, 
                              start_year    = -1,
                              end_year      = -1, 
                              compare_plots = FALSE) {
  
  # burden estimate type -- central or stochastic
  if (psa > 0) {
    burden_estimate_type <- "stochastic_burden_estimate_"
  } else {
    burden_estimate_type <- "central_burden_estimate_"
  }
  
  # cfr (case fatality rate) options (Wolfson and/or Portnoy)
  for (cfr_option in cfr_options) {
    
    # diagnostic plots filename
    pdf (paste0 (plot_folder,
                 "diagnostic_plot_",
                 antigen,
                 cfr_option, 
                 ".pdf"))
    
    # burden estimates of all scenarios
    all_burden <- NULL
    
    # scenarios
    for (scenario_name in scenarios) {
      
      # vaccine coverage file
      vaccine_coverage_file <- paste0 (vaccine_coverage_folder, 
                                       coverage_prefix, 
                                       touchstone, 
                                       antigen,
                                       scenario_name, 
                                       ".csv")
      
      # burden estimate filename
      burden_estimate_file <- paste0 (burden_estimate_type,
                                      antigen, 
                                      group_name, 
                                      scenario_name, 
                                      ".csv")
      
      
      # append/suffix cfr_option to the end of filename
      updated_burden_estimate_file <- 
        str_replace (string      = burden_estimate_file, 
                     pattern     = ".csv", 
                     replacement = paste0 ("_", cfr_option, ".csv"))
      
      # burden file with folder
      burden_file = paste0 (burden_estimate_folder,
                            cfr_option, "/", 
                            updated_burden_estimate_file)
      
      # read data -- vaccine coverage and burden estimates
      vaccine_coverage <- fread (vaccine_coverage_file)
      burden_estimate  <- fread (burden_file)
      
      # set start and end year if not set
      if (start_year == -1) {
        start_year <- min (burden_estimate [, year])
      }
      if (end_year == -1) {
        end_year <- max (burden_estimate [, year])
      }
      
      # extract vaccine coverage and burden estimates between start and end years
      vaccine_coverage <- vaccine_coverage [year >= start_year & year <= end_year]
      burden_estimate  <- burden_estimate  [year >= start_year & year <= end_year]
      
      # set vaccine (type) to SIA for campaign immunisation
      vaccine_coverage [activity_type == "campaign", vaccine := "SIA"]
      
      # add scenario name to burden estimate data table
      burden_estimate [, scenario := scenario_name]
      
      # ------------------------------------------------------------------------
      # If comparative plots across all scenarios is desired, than combine
      # burden estimates across all scenario. Make sure the combined data table
      # size is within the size of RAM.
      if (compare_plots) {
        # combine burden estimates of all scenarios
        if (is.null (all_burden)) {
          all_burden <- burden_estimate
        } else {
          all_burden <- rbindlist (list (all_burden, 
                                         burden_estimate), 
                                   use.names = TRUE)
        }
      }
      # ------------------------------------------------------------------------
      
      # if countries are specified to all, then set countries to all countries in coverage file
      if (countries == "all") {
        countries	<- as.character (unique (burden_estimate [, country] ) )  
      }
      
      # iso3 country codes
      country_iso3_codes <- countries
      
      # plot for each country
      for (country_iso3_code in country_iso3_codes) {
        
        # plot vaccine coverage
        coverage_plot <- ggplot (data = vaccine_coverage [country_code == country_iso3_code], 
                                 aes (x = year,
                                      y = coverage * 100, 
                                      color = factor (vaccine))) +  
          scale_x_continuous (breaks = pretty_breaks ()) + 
          geom_point () + 
          labs (title = countrycode (sourcevar   = country_iso3_code, 
                                     origin      = "iso3c", 
                                     destination = "country.name"),
                x = "Year", 
                y = "Vaccine coverage (%)",
                colour = "vaccine") + 
          theme_bw ()
        
        
        # plot burden -- cases, deaths, dalys
        plotwhat       <- c("cases", "deaths", "dalys")
        plotwhat_label <- c("Cases", "Deaths", "DALYs")
        
        burden_plot_list <- lapply (1:length(plotwhat), function (i) {
          
          toplot = plotwhat[i]
          
          p <- ggplot(burden_estimate [country == country_iso3_code], 
                      aes(x = year, y = get(toplot))) +
            scale_x_continuous (breaks = pretty_breaks ()) + 
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) + 
            labs (title = countrycode (sourcevar   = country_iso3_code, 
                                       origin      = "iso3c", 
                                       destination = "country.name"),
                  x = "Year", 
                  y = plotwhat_label [i]) + 
            theme_bw ()
        })
        
        # list of plots
        plot_list <- list (coverage_plot, 
                           burden_plot_list [[1]], 
                           burden_plot_list [[2]], 
                           burden_plot_list [[3]])
        
        # arrange plots in a single page
        plots <- ggarrange (plotlist = plot_list, 
                            ncol = 2, 
                            nrow = 2)
        
        # print plots
        print (annotate_figure (plots,
                                top = text_grob (scenario_name,
                                                 color = "black",
                                                 size = 9)))
        
      } # end of loop -- for (country_iso3_code in country_iso3_codes)
      
    } # end of loop -- for (scenario_name in scenarios)
    
    
    # --------------------------------------------------------------------------
    # comparative plots across all scenarios
    if (compare_plots) {
      
      # add comparative plot across all scenarios for each country
      for (country_iso3_code in country_iso3_codes) {
        
        # plot burden -- cases, deaths, dalys
        plotwhat       <- c("cases", "deaths", "dalys")
        plotwhat_label <- c("Cases", "Deaths", "DALYs")
        
        for (i in 1:length (plotwhat)) {
          
          toplot = plotwhat [i]
          
          p <- ggplot(all_burden [country == country_iso3_code], 
                      aes(x     = year, 
                          y     = get (toplot), 
                          group = scenario, 
                          color = factor (scenario))) +
            scale_x_continuous (breaks = pretty_breaks ()) + 
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) + 
            labs (title = countrycode (sourcevar   = country_iso3_code, 
                                       origin      = "iso3c", 
                                       destination = "country.name"),
                  x = "Year", 
                  y = plotwhat_label [i], 
                  colour = "Scenario") + 
            theme_bw ()
          
          print (p)
        }
      } # end of loop -- for (country_iso3_code in country_iso3_codes)
      # --------------------------------------------------------------------------
      
      
      # --------------------------------------------------------------------------
      # add comparative plot across all scenarios for each country 
      # in comparison to a base scenario, (i.e.) proportional increase or decrease
      # in cases, deaths and dalys in comparison to base scenario
      
      # burden columns
      burden_columns <- c("cases", "deaths", "dalys")
      
      # streamline burden estimates for comparison
      all_burden_compare <- all_burden [, lapply (.SD, sum, na.rm=TRUE), 
                                        .SDcols = burden_columns,
                                        by = .(country, year, scenario)]
      
      # sort by specific columns
      setorderv (all_burden_compare, c("scenario", "country", "year"))
      
      # extract burden estimates of base scenario
      burden_base <- all_burden_compare [scenario == base_scenario]
      
      all_burden_compare [, cases  := cases  / burden_base [, cases  ]]
      all_burden_compare [, deaths := deaths / burden_base [, deaths ]]
      all_burden_compare [, dalys  := dalys  / burden_base [, dalys  ]]
      
      
      for (country_iso3_code in country_iso3_codes) {
        
        # plot burden -- cases, deaths, dalys
        plotwhat       <- c("cases", "deaths", "dalys")
        plotwhat_label <- c("Cases", "Deaths", "DALYs")
        plotwhat_label <- paste0 (plotwhat_label, 
                                  " (proportion in comparison to base scenario)")
        
        for (i in 1:length (plotwhat)) {
          
          toplot = plotwhat [i]
          
          p <- ggplot(all_burden_compare [country == country_iso3_code], 
                      aes(x     = year, 
                          y     = get (toplot), 
                          group = scenario, 
                          color = factor (scenario))) +
            scale_x_continuous (breaks = pretty_breaks ()) + 
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) + 
            labs (title = countrycode (sourcevar   = country_iso3_code, 
                                       origin      = "iso3c", 
                                       destination = "country.name"),
                  x = "Year", 
                  y = plotwhat_label [i], 
                  colour = "Scenario") + 
            theme_bw ()
          
          print (p)
        }
        
        
      } # end of loop -- for (country_iso3_code in country_iso3_codes)
      # --------------------------------------------------------------------------
      
      
      # --------------------------------------------------------------------------
      # add comparative plot across all scenarios for each country 
      # total burden in each scenario for a given number of years
      # total cases, deaths and dalys for each scenario
      
      # burden columns
      burden_columns <- c("cases", "deaths", "dalys")
      
      # streamline burden estimates for comparison
      all_burden_compare <- all_burden [, lapply (.SD, sum, na.rm=TRUE), 
                                        .SDcols = burden_columns,
                                        by = .(country, scenario)]
      
      # sort by specific columns
      setorderv (all_burden_compare, c("country", "scenario"))
      
      for (country_iso3_code in country_iso3_codes) {
        
        # plot burden -- cases, deaths, dalys
        plotwhat       <- c("cases", "deaths", "dalys")
        plotwhat_label <- c("Cases", "Deaths", "DALYs")
        
        for (i in 1:length (plotwhat)) {
          
          toplot = plotwhat [i]
          
          p <- ggplot(all_burden_compare [country == country_iso3_code], 
                      aes(x     = scenario, 
                          y     = get (toplot), 
                          fill = factor (scenario))) +
            geom_col() +
            ylab (toplot) +
            labs (title = countrycode (sourcevar   = country_iso3_code, 
                                       origin      = "iso3c", 
                                       destination = "country.name"),
                  x = "Scenario", 
                  y = plotwhat_label [i], 
                  colour = "Scenario") + 
            theme_bw () +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
            expand_limits (y = 0)
          
          print (p)
        }
        
        
      } # end of loop -- for (country_iso3_code in country_iso3_codes)
      # --------------------------------------------------------------------------
      
      
      # --------------------------------------------------------------------------
      # add comparative plot across all scenarios for each country
      # in comparison to a base scenario, (i.e.) proportional increase or decrease
      # total cases, deaths and dalys for each scenario in comparison to base scenario
      
      # burden columns
      burden_columns <- c("cases", "deaths", "dalys")
      
      # streamline burden estimates for comparison
      all_burden_compare <- all_burden [, lapply (.SD, sum, na.rm=TRUE), 
                                        .SDcols = burden_columns,
                                        by = .(country, scenario)]
      
      # sort by specific columns
      setorderv (all_burden_compare, c("scenario", "country"))
      
      # extract burden estimates of base scenario
      burden_base <- all_burden_compare [scenario == base_scenario]
      
      all_burden_compare [, cases  := cases  / burden_base [, cases  ]]
      all_burden_compare [, deaths := deaths / burden_base [, deaths ]]
      all_burden_compare [, dalys  := dalys  / burden_base [, dalys  ]]
      
      for (country_iso3_code in country_iso3_codes) {
        
        # plot burden -- cases, deaths, dalys
        plotwhat       <- c("cases", "deaths", "dalys")
        plotwhat_label <- c("Cases", "Deaths", "DALYs")
        plotwhat_label <- paste0 (plotwhat_label, 
                                  " (proportion in comparison to base scenario)")
        
        for (i in 1:length (plotwhat)) {
          
          toplot = plotwhat [i]
          
          p <- ggplot(all_burden_compare [country == country_iso3_code], 
                      aes(x     = scenario, 
                          y     = get (toplot), 
                          fill = factor (scenario))) +
            geom_col() +
            ylab (toplot) +
            labs (title = countrycode (sourcevar   = country_iso3_code, 
                                       origin      = "iso3c", 
                                       destination = "country.name"),
                  x = "Scenario", 
                  y = plotwhat_label [i], 
                  colour = "Scenario") + 
            theme_bw () +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
          
          print (p)
        }
        
        
      } # end of loop -- for (country_iso3_code in country_iso3_codes)
      
    } # end of -- if (compare_plots) 
    # --------------------------------------------------------------------------
    
    dev.off ()
    
  } # end of loop -- for (cfr_option in cfr_options)
  
} # end of function -- diagnostic_plots
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Create vaccine coverage file (0% coverage) for no vaccination scenario using 
# another vaccination scenario. This is done because VIMC vaccine coverage file
# for no vaccination scenario is empty.
create_no_vaccination_coverage_file <- function (no_vaccination_coverage_file, 
                                                 vaccination_coverage_file) {
  
  # read vaccine coverage from vaccination scenario coverage file
  vaccine_coverage <- fread (file = vaccination_coverage_file)
  
  # set vaccine coverage to zero
  vaccine_coverage [, coverage := 0]
  
  # save zero vaccine coverage to no vaccination scenario coverage file
  fwrite (x    = vaccine_coverage, 
          file = no_vaccination_coverage_file)
  
  return ()
  
} # end of function -- create_no_vaccination_coverage_file
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Create vaccine coverage file for campaign only vaccination scenario using 
# (routine + campaign) vaccination scenario. 
# Set routine coverage to zero. This is done because routine coverage values are
# needed even if they are only zero to run campaign only vaccination scenario.
create_campaign_vaccination_coverage_file <- function (campaign_only_vaccination_coverage_file, 
                                                       routine_campaign_vaccination_coverage_file) {
  
  # read vaccine coverage from (routine + campaign) vaccination coverage file
  vaccine_coverage <- fread (file = routine_campaign_vaccination_coverage_file)
  
  # set vaccine coverage of routine to zero
  vaccine_coverage [activity_type == "routine", coverage := 0]
  
  # save vaccine coverage to campaign only vaccination coverage file
  fwrite (x    = vaccine_coverage, 
          file = campaign_only_vaccination_coverage_file)
  
  return ()
  
} # end of function -- create_campaign_vaccination_coverage_file
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create latin hyper cube sample of input parameters for 
# probabilistic sensitivity analysis
# ------------------------------------------------------------------------------
CreatePSA_Data <- function (psa             = 0, 
                            seed_state      = 1,
                            psadat_filename = NA) {
  
  # set seed for reproducibility
  set.seed (seed = seed_state)  
  
  # ----------------------------------------------------------------------------
  # create data table to store psa values 
  # input parameters for sensitivity analysis
  
  # run_id
  psadat <- data.table (run_id = c (1:psa))  
  
  # initialise parameters
  
  # intercept of regression function for vaccine efficacy by age
  psadat [, vaceffbyage_a := 0] 
  
  # slope of regression function for vaccine efficacy by age
  psadat [, vaceffbyage_b := 0]
  
  # vaccine efficacy (dose 2)
  psadat [, take3_input := 0]
  
  # proportional change in case fatality rate
  psadat [, mortality_input := 0]
  # ----------------------------------------------------------------------------
  
  # construct a random Latin hypercube design
  cube <- randomLHS (n = psa, 
                     # k = (ncol (psadat) - 1)
                     # vaceffbyage_a, vaceffbyage_b, and take3_input will vary
                     # together, since a, b refers to MCV1 and take3 refers to MCV2
                     k = 2  
  )
  
  # ----------------------------------------------------------------------------
  # ! vaccine efficacy parameters for function vaceffByAge (in the fortran code)
  # ve_a = real(0.64598,dp)
  # ve_b = real(0.01485,dp)
  
  # Coefficients:
  #                   Estimate  Std. Error t value Pr(>|t|)    
  # (Intercept)       0.645985   0.058727  11.000 3.49e-15 ***
  # age_at_first_dose 0.014853   0.004969   2.989  0.00426 ** 
  #
  # (Intercept)      : 0.645985 (95% CI: 0.5281395,   0.7638297)
  # age_at_first_dose: 0.014853 (95% CI: 0.004882493, 0.02482367)
  
  # intercept of regression function for vaccine efficacy by age
  mean_intercept <- 0.645985
  low_95CI       <- 0.5281395
  sd_intercept   <- (mean_intercept - low_95CI) / 1.96
  
  psadat [, vaceffbyage_a := qtruncnorm (cube [, 1], 
                                         a    = (mean_intercept - 3 * sd_intercept), 
                                         b    = (mean_intercept + 3 * sd_intercept), 
                                         mean = mean_intercept, 
                                         sd   = sd_intercept 
  ) ]
  
  # slope of regression function for vaccine efficacy by age
  mean_slope <- 0.014853
  low_95CI   <- 0.004882493
  sd_slope   <- (mean_slope - low_95CI) / 1.96
  
  psadat [, vaceffbyage_b := qtruncnorm (cube [, 1], 
                                         # a    = (mean_slope - 3 * sd_slope), == -0.0004079801
                                         a    = 0,
                                         b    = (mean_slope + 3 * sd_slope), 
                                         mean = mean_slope, 
                                         sd   = sd_slope
  ) ]
  # ----------------------------------------------------------------------------
  
  # vaccine efficacy (dose 2)
  # mean 98% -- +- ~ 2%
  psadat [, take3_input := qtruncnorm (cube [, 1], 
                                       a    = (0.98 - 0.02), 
                                       b    = (0.98 + 0.02), 
                                       mean = 0.98, 
                                       sd   = (0.02/3) 
  ) ]
  
  
  # proportional change in case fatality rate
  # truncated lognormal distribution for up to 25% change
  psadat [, mortality_input := qtruncnorm (cube [, 2], 
                                           a    = (1 - 0.25), 
                                           b    = (1 + 0.25), 
                                           mean = 1, 
                                           sd   = (0.25/3) 
  ) ]
  
  # rename columns of vaceffbyage_a to take1_input and vaceffbyage_b to take2_input
  # this is done take1_input and take2_input were used earlier in the fortran model
  # which meant: take1_input # vaccine efficacy (dose 1, before age 1)
  #              take2_input # vaccine efficacy (dose 1, after age 1)
  # In this psa data table and file: take1 refers to vaceffbyage_a (intercept)
  #                                  take1 refers to vaceffbyage_b (slope)
  setnames(psadat, 
           old = c("vaceffbyage_a", "vaceffbyage_b"), 
           new = c("take1_input",   "take2_input")
  )
  
  
  # save psa data file, if file name is provided
  if (!is.na (psadat_filename)) {
    fwrite (x    = psadat, 
            file = psadat_filename)
  }
  
  # return psa data table
  return (psadat)
  
} # end of function -- CreatePSA_Data
# ------------------------------------------------------------------------------



