# functions.R
# Main functions for running the DynaMICE model

# ------------------------------------------------------------------------------
# Measles vaccine impact model
#   To estimate the health impact of measles vaccination for a given set of
#   countries.
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Execute the fortran model for a country and a PSA run
#'
#' A function nested under \code{\link{runScenario}} to run the fortran codes
#'  for measles vaccination, given a particular country and a variable set of
#'  probabilistic sensitivity analysis (PSA).
# ------------------------------------------------------------------------------
# Changes 2019:
# 1) expanded matrix is rescaled to keep the total number of
# contacts in weekly age groups correct (checked - R0 of rescaled matrix is now
# the same as prior to expanding it to weekly ages)
# 2) projects the contact matrix to represent country's demography
# ------------------------------------------------------------------------------
#' @param ii Numeric order of the selected country for evaluation.
#' @param iso3 ISO-3 code of the selected country.
#' @param years A vector containing continuous calender years for simulation.
#' @param vaccination A Numeric indicator that determines vaccination programmes
#'  for children: 0 - No vaccination, 1 - Only MCV1,  2 - MCV1 and MCV2.
#' @param using_sia A numeric indicator that determines Whether supplementary
#' immunisation activities (SIAs) are implemented: 0 - no SIA, 1 - with SIA.
#' @param dinf Duration of infection in days.
#' @param gamma_rate Recovery rate, with a unit of 1/(\code{dinf}*\code{tstep}*year).
#' @param r0_basic Basic reproduction number (R0) calculated by the eigenvalue
#' of the product of contact matrix (\code{contact}) and duration of infection
#' (\code{dinf}).
#' @param amplitude A numeric variable for amplitude for seasonality.
#' @param take A 1*3 vector that denotes the level of vaccine protection using
#' the 'take' ('all-or-not') assumption. Vaccine protection is assumed to vary
#' by age under a linear regression, with (1) and (2) elements of the vector
#' representing the intercept and coefficient of the first dose, respectively.
#' (3) element represents the vaccine protection by the second dose. Note that
#' in the previous fortran model without the incorporation of age-related
#' regression model, (1) and (2) elements were directly indicating protection of
#'  the first dose before and after 1 year old, respectively.
#' @param degree A 1*3 vector that denotes the level of vaccine protection using
#'  the 'degree' ('leak') assumption. (1) and (2) elements represents the
#'  protection of the first dose before and after 1 year old, respectively, and
#' (3) element shows that of the second dose. Note that the second dose only has
#'  an effect if \code{vaccine=2}.
#' @param sia.method A Numeric indicator that determines how SIAs are carried
#' out: 1 - varying \code{take}, 2 - varying \code{degree}.
#' @param coverage_routine A data frame for routine vaccination coverage under
#' a selected scenario.
#' @param coverage_sia A data frame for SIA coverage under a selected scenario.
#' @param timeliness A data frame for timeliness estimates by age.
#' @param contact A data frame for the contact matrix by age (0-100 years old).
#' @param c_rnought A numeric variable of R0 value for the selected country
#' \code{iso3}. Estimates provided by Ferrari and colleagues.
#' @param population A data frame for population size by country and age.
#' @param tstep Simulation time steps for each year.
#' @param save.scenario A folder name of the selected scenario.
#' @param foldername A folder name for input and output files of the selected
#' scenario. It typically does not exist but may come in handy when only
#' processing results.
#' @param measles_model An executable file that processes fortan codes of the
#' measles model.
#' @param debug_model A logical variable that determines whether to debug the
#' model.
#' @param debug_spinup A logical variable that determines whether to generate
#' outputs for the spin-up (equilibrium) period.
#' @param debug_age A numeric variable that determines the format of output age
#' groups: 0 - all in annual age-strata, 1 - age between 0 and 2 in weekly
#' age-strata and age between 3 and 100 in annual age-strata.
#' @param debug_compartment A numeric variable that determines the compartments
#' to output: 0 - number of cases, and 1 - each compartment.
#' @param debug_country A vector including ISO3 codes of country to debug, or *
#' to debug all countries.
#' @param debug_relative A logical variable that determines the output format of
#' new cases: TRUE - proportion, and FALSE - absolute number.
#' @param debug_timepoints A logical variable that determines the output unit of
#' time: 0 - per year, 1 - per timepoint and report first 25% of timepoints, and
#' 2 - per timepoint and report all timepoints.
#' @param r A numeric variable of the order of the PSA runs.
#' @param runs A numeric variable of the total runs for each scenario.
#' @param psa A numeric variable of the total runs for PSA. Use 0 to indicate a
#' single run without PSA.
#' @param psa_var A data frame including the parameter values used in PSA.
#' @param run_model A logical variable that determines whether to run the
#' measles model.
#' @param remove_files A logical variable that determines whether to remove
#' output files after finishing a single country run.
#' @param contact_mat A character variable that indicates the assumption for
#' contact pattern. "nomix" - homogeneous mixing, "uk" - POLYMOD Great Britain
#' physical contacts, and "syn" - country-specific synthetic matrix.
#' @examples
#' runCountry (
#'   ii                 = 2,
#'   iso3               = "BGD",
#'   years              = 1980:2020,
#'   vaccination        = 1,
#'   using_sia          = 1,
#'   dinf               = 14,
#'   gamma_rate         = 0.02607,
#'   r0_basic           = 13,
#'   amplitude          = 0.05,
#'   take               = c(0.64598, 0.01485, 0.98),
#'   degree             = c(0.85, 0.95, 0.98),
#'   sia.method         = 1,
#'   coverage_routine   = coverage_routine,
#'   coverage_sia       = coverage_sia,
#'   timeliness         = timeliness,
#'   contact            = contact,
#'   c_rnought          = 10,
#'   population         = population,
#'   tstep              = 1000,
#'   save.scenario      = "scenario01",
#'   foldername         = NULL,
#'   measles_model      = "vaccine2019_sia_singlematrix.exe",
#'   debug_model        = FALSE,
#'   debug_spinup       = FALSE,
#'   debug_age          = 0,
#'   debug_compartments = 1,
#'   debug_country      = "*",
#'   debug_relative     = FALSE,
#'   debug_timepoints   = 0,
#'   r                  = 1,
#'   runs               = 1,
#'   psa                = 0,
#'   psa_var            = psa_var,
#'   run_model          = TRUE,
#'   remove_files       = FALSE,
#'   contact_mat        = "nomix")
runCountry <- function (
  #variables specific for loop
  ii,
  iso3,
  years,

  #infection dynamic variables
  vaccination,
  using_sia,
  dinf,
  gamma_rate,
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
  c_rnought,
  population,

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
  run_model,
  remove_files,
  contact_mat
) {

  # temporarily assign 3-letter-ISO code to Kosovo until Kosovo is assigned official ISO3-code
  if (iso3 == "XK"){
    fortran_country_code <- "XKX"
  } else {
    fortran_country_code <- iso3
  }

  if (psa > 0) {
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

    # remove 0% coverage SIAs
    coverage_sia <- setDT (coverage_sia) [country_code == iso3 & coverage != 0]

    if ((using_sia == 1) & (nrow (coverage_sia) > 0)) {

      sia_a0       <- coverage_sia [, a0]
      sia_year     <- coverage_sia [, year]
      sia_a1       <- coverage_sia [, a1]
      sia_coverage <- coverage_sia [, coverage]

      sia_coverage [which (sia_coverage>1)] <- 0.95

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

    # get timepoints for SIA years, assuming that year 1 of interest is the year 1980
    # assume SIAs to take place at the beginning of a calendar year,
    # with 1-month interval (1000/12 timesteps) for multiple SIA rounds within a single year

    t_sia <- unlist (sapply (sort (unique (sia_year)), function (iyr, sia_year)
      return (t_run[1] + tstep * (iyr - years[1]) + c(1:sum(sia_year == iyr) - 1)*(round(1000/12))),
      sia_year = sia_year
      ))

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

    # country specific timeliness curve
    country_timeliness <- timeliness [country_code == iso3 & !is.na(age), timeliness]
    timeliness_ages    <- timeliness [country_code == iso3 & !is.na(age), age]

    # country_specific contact matrix
    q             <- c_rnought / r0_basic  # proportionality factor (infectivity, underreporting)
    contact_day   <- contact * q
    contact_tstep <- contact_day * (365 / tstep)
    r0_tstep  <- Re (eigen (contact_tstep, only.values=T)$values[1]) # will need it to make sure R0 is the same after rescaling

    # Beta only a single file
    s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
    jt        <- 3  # how many ages to expand to s (or to weeks)
    beta_full <- matrix (0,
                        ncol = 254,
                        nrow = 254)

    # expand contact of first 3 age-strata with itself
    beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
      A = contact_tstep [1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same
      expand_rows =  s, expand_cols =  s,
      rescale_rows = FALSE, rescale_cols = FALSE)

    # expand contact for first 3 age-strata with all other contacts
    beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
      A = contact_tstep [1:jt,(jt+1):ncol(contact_tstep)],
      expand_rows = s, expand_cols = 1,
      rescale_rows = F, rescale_cols = F)

    beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
      A = contact_tstep [(jt+1):nrow(contact_tstep),1:jt]/s,  # adjust to ensure the mean total number of contacts stays the same
      expand_rows = 1, expand_cols = s,
      rescale_rows = F, rescale_cols = F)

    beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
      contact_tstep [(jt+1):nrow(contact_tstep),(jt+1):ncol(contact_tstep)]

    # make the matrix reciprocal
    beta_full <- (beta_full + t(beta_full)) / 2
    # [method used previously] adjust using age distribution of population
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

      #old script groups those aged 70-80, but division is by actual population size
      pop.vector <- population[country_code == iso3 & year == y, value]

      # first expand polymod matrix (contact_tstep) and population vector and
      # then divide by population sizes, otherwise it doesn't work.
      pop.vector_full <- c(rep(pop.vector[1:jt]/s, each = s), pop.vector[(jt +1):length(pop.vector)])

      # change zero values to 1 to avoid division by 0
      pop.vector_full[pop.vector_full==0] <- 1


      if (vaccination >= 1) {

        # Maximum coverage can (obviously) only be 100%
        # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
        # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
        # In essence, this becomes the inverse of the cumulative timeliness curve
        cycov <- coverage_routine [country_code == iso3 & year == y & vaccine == "MCV1", coverage] /
          timeliness [country_code == iso3 & is.na(age), prop_final_cov]

        if (length(cycov) == 0) {   # check if vaccine is not yet introduced and thereby, coverage value missing for this year
          cycov <- 0
        } else if (is.na(cycov)) {
          cycov <- 0
        }

        country_year_timeliness_mcv1 <- 1 - min(cycov,1) * country_timeliness

        country_year_timeliness_mcv1 <- -diff(country_year_timeliness_mcv1) /
          (country_year_timeliness_mcv1 [1:(length(country_year_timeliness_mcv1)-1)])

        # Timeliness is reported by week in the first year, and by month in the second year. Assume there is no vaccination in between
        country_year_timeliness_mcv1 [is.na(country_year_timeliness_mcv1)] <- 0
        country_year_timeliness_mcv1 [is.nan(country_year_timeliness_mcv1)] <- 0
        country_year_timeliness_mcv1_allages <- rep(0, 254)
        country_year_timeliness_mcv1_allages [round(timeliness_ages)] <- country_year_timeliness_mcv1

      } else {
        country_year_timeliness_mcv1_allages <- rep(0, 254)
      }

      if(vaccination == 2){
        country_year_mcv2 <- coverage_routine [country_code == iso3 & year == y & vaccine == "MCV2", coverage]
      } else {
        country_year_mcv2 <- 0
      }

      if ( length (country_year_mcv2) == 0 | is.na(country_year_mcv2) ) {
        country_year_mcv2 <- 0
      }


      # write data for each year
      iyr <- which (years == y)
      beta_path <- paste0 ("outcome/", save.scenario, "/", foldername,
                           "/input/", psadir, fortran_country_code)

      if (iyr<10) {
        dynamice_input_file <- paste0 (beta_path, "_00", iyr, "measle_data.txt")

      } else if(iyr<100) {
        dynamice_input_file <- paste0 (beta_path,  "_0", iyr, "measle_data.txt")

      } else {
        dynamice_input_file <- paste0 (beta_path,   "_", iyr, "measle_data.txt")
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
        gamma_rate,
        tstep,
        take1,
        take2,
        take3,
        degree1,
        degree2,
        degree3,
        length(sia_a0),
        sia_a0,
        c_rnought,
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
          pop.vector_full,
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
} # end of function -- runCountry
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Run the measles model for a selected vaccination scenario
#'
#' A function that execute the measles model under a selected vaccination
#' scenario, including a pre-specified set of countries and runs of probability
#' sensitivity analysis (PSA).
#'
#' @param vaccine_coverage_folder A folder name for the vaccine coverage files.
#' @param coverage_prefix A prefix used in the name of vaccine coverage file.
#' @param touchstone A version note in the file name used by VIMC. Include a
#' underscore at the beginning and end of the name.
#' @param antigen A disease name used by VIMC: "Measles".
#' @param scenario_name A name of vaccination scenarios.
#' @param scenario_number A folder name of the selected scenario.
#' @param vaccine_coverage_subfolder A folder name under the \code{x} folder for
#' the vaccine coverage files.
#' @param burden_estimate_folder A folder name for the file which contains the
#' model outputs for evaluation. Include a slash at the end.
#' @param group_name A modelling group name used by VIMC.
#' @param countries A vector of ISO-3 country codes used in the analysis. Use
#' "all" to include all countries.
#' @param cluster_cores A number of cores to be used in the cluster.
#' @param psa A numeric variable of the total runs for PSA. Use 0 to indicate a
#' single run without PSA.
#' @param vaccination A numeric indicator that determines vaccination programmes
#'  for children: 0 - No vaccination, 1 - Only MCV1,  2 - MCV1 and MCV2.
#' @param using_sia A numeric indicator that determines Whether supplementary
#' immunisation activities (SIAs) are implemented: 0 - no SIA, 1 - with SIA.
#' @param measles_model An executable file that processes Fortran codes of the
#' measles model.
#' @param debug_model A logical variable that determines whether to debug the
#' model.
#' @param contact_mat A character variable that indicates the assumption for
#' contact pattern. "nomix" - homogeneous mixing, "uk" - POLYMOD Great Britain
#' physical contacts, and "syn" - country-specific synthetic matrix.
#' @param step_ve A logical variable that determines whether a step change on
#' age-related vaccine efficacy is applied. This implies different meaning of
#' vaccine efficacy vector (\code{take}) and a corresponding \code{measles_model}
#' should be applied.
#' @examples
#' runScenario (
#'   vaccine_coverage_folder    = "vaccine_coverage_upd/",
#'   coverage_prefix            = "coverage",
#'   touchstone                 = "_201910gavi-5_",
#'   antigen                    = "measles-",
#'   scenario_name              = "campaign-only-default",
#'   scenario_number            = scenario_number,
#'   vaccine_coverage_subfolder = "scenarios/"
#'   burden_estimate_folder     = "central_burden_estimate/",
#'   group_name                 = "LSHTM-Jit-",
#'   countries                  = c("BGD","ETH"),
#'   cluster_cores              = 3,
#'   psa                        = 0,
#'   vaccination                = 0,
#'   using_sia                  = 1,
#'   measles_model              = "vaccine2019_sia_singlematrix.exe",
#'   debug_model                = FALSE,
#'   contact_mat                = "nomix",
#'   step_ve                    = FALSE)
runScenario <- function (vaccine_coverage_folder    = "",
                         coverage_prefix            = "",
                         touchstone                 = "",
                         antigen                    = "",
                         scenario_name,
                         scenario_number,
                         vaccine_coverage_subfolder = "",
                         burden_estimate_folder,              # burden estimate folder
                         group_name,                          # modelling group name
                         countries                  = "all",
                         cluster_cores              = 1,
                         psa                        = 0,      # psa runs; 0 for single run
                         vaccination,  # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
                         using_sia,    # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA
                         measles_model,                       # measles model
                         debug_model                = FALSE,  # debug model (T/F)
                         contact_mat,                         # contact matrix: "nomix","uk","syn"
                         step_ve                              # step change of take
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
  # input data file names
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

  # data file name for PSA variables
  # same for each scenario. should be CREATED before simulation using CreatePSA_Data()
  data_psa 		<- "psa_variables.csv"

  # --------------------------------------------------------------------------


  # --------------------------------------------------------------------------
  # Advanced options
  # --------------------------------------------------------------------------

  sia.method	<- 1				             # 1 = variable take; 2 = variable degree
  dinf		    <- 14				             # duration of infection (days)
  amplitude 	<- 0.05			             # amplitude for seasonality

  # tack [3] : vaccine efficacy (dose2) and only has an effect if vaccine==2.
  # take [1:2] have different definitions, and should be used with corresponding measles models.

  if (step_ve) {

    take 		  <- c (0.85, 0.95, 0.98)
    # take [1] # vaccine efficacy (dose 1, before age 1)
    # take [2] # vaccine efficacy (dose 1, after age 1)

  } else {

    take 		  <- c (0.64598, 0.01485, 0.98)
    # take [1] refers to vaceffbyage_a (intercept)
    # take [2] refers to vaceffbyage_b (slope)
  }

  degree 		  <- c (0.85, 0.95, 0.98)    # vaccine efficacy for degree1 (degree dose 1, before age 1), degree2 (dose 1, after age 1) & degree3 (dose 2). Note that dose2 only has an effect if vaccine==2.
  tstep			  <- 1000				             # Number of time steps in a year
  gamma_rate  <- 1 / (dinf * tstep/365)  # rate of losing infectivity

  # ----------------------------------------------------------------------------
  # Measles model
  # filename of compiled fortran-model (in ./model/compiled/)
  # ----------------------------------------------------------------------------
  # measles_model <- "vaccine2019_sia_singlematrix.exe" # change this to reflect the right version of the fortran code
  measles_model <- measles_model
  # ----------------------------------------------------------------------------


  # number of clusters to use
  # if larger than 1, country-specific model runs are distributed over specified number of clusters
  # note that model uses a lot of memory, so might not want to max out all clusters
  use_cluster  <- cluster_cores   # debug #
  remove_files <- TRUE
  run_model    <- TRUE
  # folder will be created if not given - should usually be commented out, except when run model is FALSE

  # ------------------------------------------------------------------------------
  # Debug
  # ------------------------------------------------------------------------------

  debug_country		  <- "*"		    	#ISO3 codes of country to debug, * to debug all countries
  debug_spinup		  <- FALSE	    	#TRUE/FALSE: If true, generate data for spin-up period of model
  debug_model       <- debug_model	#TRUE/FALSE: If true: generate data for period after spin-up
  debug_compartments<- 1    	      #0-2: If 0: output number of cases. If 1: output size of each compartment. If 2: debug vaccinated
  debug_age         <- 0            #0-2. If 0: output all in annual age-strata. If 1: output age 0-2 in weekly age-strata, 3-100 in annual age-strata. If 2: sum all age-strata.
  debug_timepoints	<- 0			      #0-2. If 0: output per year. If 1: output per timepoint and report first 25% of timepoints. If 2: output per timepoint and report all timepoints.
  debug_relative		<- FALSE	    	#If true: output proportion of new cases. If false, output absolute number of new cases.


  # START OF MODEL
  # ----------------------------------------------------------------------------
  # SETUP

  # load correct libraries when on cluster (using open MPI)
  if ("cluster_using_openmpi" %in% commandArgs()){
    using_openmpi <- TRUE
    require(doMPI)
    require(parallel)
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
  timeliness  		  <- setDT (data_timeliness)
  rnought	    		  <- setDT (data_r0)
  population  		  <- setDT (data_pop)
  template    		  <- setDT (data_template)

  # --------------------------------------------------------------------------
  # if psa variables file does not exist, then create psa variables file
  if (psa > 0) {
    if (file.exists (data_psa)) {

      # read csv if file already exists
      psa_var <- fread (data_psa)

      # check if psa_var corresponds with psa
      if(nrow(psa_var) != psa){
        stop(paste0("Number of runs in PSA file is not the same as those specified! Variables used in the same run in each scenario should be similar. Delete or rename ",data_psa," if a new file needs to be created, or change the number of PSA runs to ",nrow(psa_var)," if the same file should be used."))
      }
    } else {
      stop(paste0("There is no files for inputting PSA variables. Use CreatePSA_Data to generate the file."))
      }
  }
  # --------------------------------------------------------------------------

  # process data

  # if countries are specified to all, then set countries to all countries in coverage file
  if (countries[[1]] == "all") {
    countries	<- as.character (unique (coverage_routine [, country_code] ) )
  }

  # start and end years (should go in as input to function -- INPUT-FUNCTION)
  years <- as.numeric (c(1980:2100))

  for (cty in countries) {
    if (psa > 1){
      for (p in 1:psa){
        write(
          paste0(
            cty,
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
          cty,
          " ",
          "1 0"
        ),
        file   = "gavi_progress",
        append = TRUE
      )
    }
  }

  # Process contact matrices
  if (contact_mat == "syn"){
    contact_list  <- sapply (countries,
                             function(cty){data_contact_syn[[cty]]},
                             simplify = FALSE, USE.NAMES = TRUE)
    r0_basic_list <- sapply (countries,
                             function(cty){Re (eigen (contact_list[[cty]] * dinf,
                                                      only.values = T)$values[1])},
                             simplify = TRUE, USE.NAMES = TRUE)
    }

  if (contact_mat == "uk"){
    ctmat         <- as.matrix (data_contact_uk)
    r0_basic      <- Re (eigen (ctmat * dinf, only.values = T)$values[1])
    contact_list  <- sapply (countries,
                             function(x = NULL){copy(ctmat)},
                             simplify = FALSE, USE.NAMES = TRUE)
    r0_basic_list <- sapply (countries,
                             function(x = NULL){copy(r0_basic)},
                             simplify = TRUE, USE.NAMES = TRUE)
    }

  if (contact_mat == "nomix"){
    ctmat         <- matrix (1/101, nrow = 101, ncol = 101)  # eigenvalue = 1
    contact_list  <- sapply (countries,
                             function(x = NULL){copy(ctmat)},
                             simplify = FALSE, USE.NAMES = TRUE)
    r0_basic_list <- sapply (countries,
                             function(x = NULL){dinf},        # R0 = beta*dinf
                             simplify = TRUE, USE.NAMES = TRUE)
    }


  # Run model
  writelog ("gavi_log", paste0 ("Main; Foldername: ", foldername))

  if (using_openmpi | use_cluster > 1){
    if (!using_openmpi){

      if (use_cluster != parallel::detectCores()) {
        warning (paste0(parallel::detectCores(), " cores detected but ", use_cluster, " specified."))
      }

      cl <- parallel::makeCluster (use_cluster)
      doParallel::registerDoParallel (cl)

    } else {
      print (paste0 ("Clustersize: ", doMPI::clusterSize(cl)))
    }
  }


  # foreach will run countries and PSA runs in parallel if a parallel backend
  # is registered, and sequentially otherwise
  require(foreach)
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
                           iso3               = countries[ii],
                           years              = years,
                           vaccination        = vaccination,
                           using_sia          = using_sia,
                           dinf               = dinf,
                           gamma_rate         = gamma_rate,
                           r0_basic           = r0_basic_list[countries[ii]],
                           amplitude          = amplitude,
                           take               = take,
                           degree             = degree,
                           sia.method         = sia.method,
                           coverage_routine   = coverage_routine,
                           coverage_sia       = coverage_sia,
                           timeliness         = timeliness,
                           contact            = contact_list[[countries[ii]]],
                           c_rnought          = rnought[country_code == countries[ii], r0],
                           population         = population,
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
                           run_model          = run_model,
                           remove_files       = remove_files,
                           contact_mat        = contact_mat
    )
    return(out_run)
  }

  if(using_openmpi | use_cluster > 1){
    if (!using_openmpi){
      parallel::stopCluster(cl)
    } else {
      doMPI::closeCluster(cl)
    }
  }

  # check for errors
  errorcount <- 0

  for (i in 1:length(combine)){
    if ("error" %in% class(combine[[i]])){
      errormessage <- paste0("Error in task ",i,": ",combine[[i]])
      warning(errormessage)
      writelog (paste0 (gavi.dir,"gavi_log"),errormessage)
      errorcount <- errorcount + 1
      #remove from data
      combine[[i]] <- NULL
    }
  }
  # if(errorcount > 0){
  # stop(paste0("There were ",errorcount," errors."))
  # }
  # ----------------------------------------------------------------------------


  # process results

  report_years <- sort (unique (template$year))
  ages         <- sort (unique (template$age))

  # get years that the model was run for (can be different from reporting years)
  years 		   <- sort (unique (as.numeric (coverage_routine[,year])))


  # ANALYSE RUNS

  # get the names of all the files in all the subfolders the stochastic runs
  # specified in scenarioname; separate cases and popsize files
  # will use those to read them all in at once and put in a data.table
  myfiles <- list.files (path = paste0 ("outcome/", save.scenario, "/",foldername,"/output/"),
                         recursive = T, pattern = "cases", full.names = T)
  myfiles.popsize <- list.files (path = paste0 ("outcome/", save.scenario, "/",foldername,"/output"),
                                 recursive = T, pattern = "popsize", full.names = T)

  # read in all those files specified in myfiles and put them in a single data table
  all_cases  <- rbindlist (lapply (myfiles, function (fn, ...) {
    res <- withCallingHandlers (
      fread (fn, stringsAsFactors = F, check.names = F, fill = T,
             col.names = as.character (c(0:100))),
      warning = function(w) {warning(w, fn); }
      )
    res[, country := gsub ("^.+/(\\w+)_age.+$","\\1", fn)]             # get the country code from the filename
    # get run_id for psa runs
    if (psa > 0) {
      res[, run_id := as.integer (gsub("^.+run(\\d{3}).+$","\\1", fn))]  # get the run_id from the filename (subfolder)
    }
    res[, year := years]                                               # add year of simulation (from coverage data file)
  }))

  # switch back Kosovo to XK otherwise montagu returns error on upload
  all_cases[country == "XKX", country := "XK"]

  all_popsize  <- rbindlist (lapply (myfiles.popsize, function (fn, ...) {
    res <- fread (fn, stringsAsFactors = F, check.names = F, col.names = as.character(c(0:100)))
    res[, country := gsub("^.+/(\\w+)_age.+$","\\1", fn) ]             # get the country code from the filename
    if (psa > 0) {
      res[, run_id := as.integer (gsub("^.+run(\\d{3}).+$","\\1", fn)) ]  # get the run_id from the filename (subfolder)
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
  all_cases.m <- melt (all_cases, id.vars = column_names,
                       measure.vars = c(as.character(0:100)),
                       variable.name = "age", value.name = "cases", variable.factor = F)

  # all_popsize.m <- melt(all_popsize, id.vars = c("year", "country"),
  all_popsize.m <- melt (all_popsize, id.vars = column_names,
                         measure.vars = c(as.character(0:100)),
                         variable.name = "age", value.name = "cohort_size", variable.factor = F)

  # melt produces character variable when variable.factor is set to FALSE - change it to integer
  all_cases.m   [, age := lapply(.SD, as.integer), .SDcols = "age"]
  all_popsize.m [, age := lapply(.SD, as.integer), .SDcols = "age"]

  # merge cases and cohort sizes
  all_runs = merge (all_cases.m, all_popsize.m, by = c(column_names, "age"), all.x = T)

  # add country_name, life expectancy and disease (Measles) to match template file
  country_names <- unique (subset (template, select = c("country", "country_name")))
  c_names <- country_names$country_name; names(c_names) = country_names$country

  all_runs[, c("country_name", "disease") := list (c_names[country], "Measles")]


  # ----------------------------------------------------------------------------
  # add data column for remaining life expectancy
  lexp_remain <- tailor_data_lexp_remain(countries)

  if (psa > 0) {
    all_runs <- lexp_remain [all_runs,
                                  .(run_id, i.country, year, age, cases, cohort_size, country_name, disease, value),
                                  on = .(country_code = country,
                                         age_from    <= age,
                                         age_to      >= age,
                                         year         = year) ]
  } else {
    all_runs <- lexp_remain [all_runs,
                                  .(i.country, year, age, cases, cohort_size, country_name, disease, value),
                                  on = .(country_code = country,
                                         age          = age,
                                         year         = year) ]
  }

  # rename column names for output
  setnames (all_runs,
            old = c("i.country", "value"      ),
            new = c("country"  , "remain_lexp"))

  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # MCV1 coverage
  coverage_routine_MCV1 <- coverage_routine [vaccine == "MCV1"]

  # add MCV1 column
  if (psa > 0) {
    all_runs <- coverage_routine_MCV1 [all_runs,
                                       .(run_id, i.country, i.year, age, cases, cohort_size, country_name, disease, coverage, remain_lexp),
                                       on = .(country_code = country,
                                              year         = year) ]
  } else {
    all_runs <- coverage_routine_MCV1 [all_runs,
                                       .(i.country, i.year, age, cases, cohort_size, country_name, disease,  coverage, remain_lexp),
                                       on = .(country_code = country,
                                              year         = year) ]
  }

  # rename column "coverage" to "MCV1"
  setnames (x = all_runs,
            old = c("i.country", "i.year", "coverage"),
            new = c("country"  , "year"  , "MCV1"    ))
  # ----------------------------------------------------------------------------

  # # Not in use 2021:
  # # Previous methods for calculating deaths and dalys.
  # # Additional input data for CFR and life expectancy at birth are needed to recover.
  # all_runs[, deaths := cases * cfr.value]
  # all_runs[, dalys := (cases - deaths)*0.002 + deaths*(LE0 - 4)]       # assumes all deaths are age 4
  # all_runs[, dalys := (cases - deaths) * 0.002 + deaths * (LE0 - age)] # to account for actual age of death
  # ----------------------------------------------------------------------------
  # YLL calculation above is fine since most measles burden < 5 years
  # remaining life expectancy estimation can be improved by using remaining life
  # years by age, year and country
  # ----------------------------------------------------------------------------

  # OUTPUT RUNS
  # keep all columns
  output_runs <- subset (all_runs, year %in% report_years)

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
  fwrite (x    = output_runs [order (country, year, age)],
          file = paste0 (burden_estimate_folder,  burden_estimate_file))

  # clean environment
  writelog ("gavi_log", paste0 ("Main; gavi.r finished"))

  if (using_openmpi) {
    Rmpi::mpi.quit()
  }

  # return burden estimate filename (cases)
  return (burden_estimate_file)

} # end of function -- runScenario
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Estimate deaths and DALYs of measles
#'
#' A function that estimates the number of deaths by applying case fatality
#' rates (CFRs) to the number of cases, and then calculates the
#' disability-adjusted years (DALYs)
#'
#' @param cfr_options Methods for CFR estimates: "Wolfson" or/and "Portnoy".
#' @param burden_estimate_file A file name for model outputs to be saved.
#' @param burden_estimate_folder A folder name for the file which contains the
#' model outputs for evaluation. Include a slash at the end.
#' @param vimc_scenario A scenario name for vaccine coverage used in
#' \code{data_cfr_portnoy.rda}, for extracting CFR estiamtes by the Portnoy
#' method.
#' @param portnoy_scenario A scenario name for future CFR trends when using the
#' Portnoy method: "s4" - declining after 2018, or "s6" - staying stagnant at
#' the 2018 level.
#' @param psa A numeric variable of the total runs for PSA. Use 0 to indicate a
#' single run without PSA.
#' @examples
#' estimate_deaths_dalys (
#'  cfr_option             = "Portnoy",
#'  burden_estimate_file   = burden_estimate_file,
#'  burden_estimate_folder = "central_burden_estimate/",
#'  vimc_scenario          = "campaign-only",
#'  portnoy_scenario       = "s6",
#'  psa                    = 0)
# ------------------------------------------------------------------------------
#   save results in corresponding cfr_option subfolder
#   append cfr_option to burden estimates file
# ------------------------------------------------------------------------------
estimate_deaths_dalys <- function (cfr_option,
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

    # read CFRs  -- Wolfson
    cfr = setDT (data_cfr_wolfson)

    # cfr estimates for ages below 5 years, not varying by time
    # cfrs for ages above 5 years are half the cfr values for ages below 5 years

    # add CFR data column to burden estimates
    if (psa > 0) {
      burden <- cfr [burden,
                     .(run_id, disease, year, age, i.country, country_name, cohort_size, cases, CFR, remain_lexp),
                     on = .(country_code = country) ]
    } else {
      burden <- cfr [burden,
                     .(disease, year, age, i.country, country_name, cohort_size, cases, CFR, remain_lexp),
                     on = .(country_code = country) ]
    }

    # estimate deaths
    # cfrs for ages above 5 years are half the cfr values for ages below 5 years
    burden [age <  5, deaths := cases * CFR  ]
    burden [age >= 5, deaths := cases * CFR/2]

    # rename country column
    setnames (burden, old = "i.country", new = "country")

  }
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # use CFRs -- Portnoy

  if (cfr_option == "Portnoy") {

    # read CFRs (rates between 0 and 1) -- Portnoy
    cfr = setDT (data_cfr_portnoy)

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
              old = c(paste0 ("cfr_under5_", vimc_scenario, "_", portnoy_scenario),
                       paste0 ("cfr_over5_" , vimc_scenario, "_", portnoy_scenario)),
              new = c("under5", "over5"))

    # add CFR data column to burden estimates
    if (psa > 0) {
      burden <- cfr [burden,
                     .(run_id, disease, year, age, country, country_name, cohort_size, cases, over5, under5, remain_lexp),
                     on = .(country_code = country,
                            year         = year) ]
    } else {
      burden <- cfr [burden,
                     .(disease, year, age, country, country_name, cohort_size, cases, over5, under5, remain_lexp),
                     on = .(country_code = country,
                            year         = year) ]
    }

    # estimate deaths for ages under 5 years
    burden [age < 5, deaths := cases * under5]

    # estimate deaths for ages over 5 years
    burden [age >= 5, deaths := cases * over5]
  }
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # DALYs
  # calculate dalys = (ylds) + (ylls)
  burden [, dalys := ((cases - deaths) * 0.002) + (deaths * remain_lexp)]

  # adjust columns for output
  if (psa > 0) {
    save.cols <- c("run_id", names(data_template))
  } else {
    save.cols <- names(data_template)
  }
  burden <- subset (burden, select = save.cols)
  # ----------------------------------------------------------------------------

  # append/suffix cfr_option to the end of filename
  updated_burden_estimate_file <- stringr::str_replace (string      = burden_estimate_file,
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

} # end of function -- estimate_deaths_dalys
# ------------------------------------------------------------------------------

