# utils.R
# Functions for supporting utilities in the DynaMICE model (Dynamic Measles Immunisation Calculation Engine)

# ------------------------------------------------------------------------------
#' Create vaccine coverage files by scenarios
#'
#' Generates files for routine programmes (MCV1, MCV2) and supplementary
#' immunisation activities (SIAs) from the vaccine coverage files.
#'
#' @param vaccine_coverage_folder A folder name for the vaccine coverage files.
#' Include a slash at the end.
#' @param vaccine_coverage_subfolder A folder name under the \code{x} folder for
#' the vaccine coverage files.
#' @param coverage_prefix A prefix used in the name of vaccine coverage file.
#' @param touchstone A version note in the file name used by VIMC. Include a
#' underscore at the beginning and end of the name.
#' @param antigen A disease name used by VIMC: "Measles".
#' @param scenario_name A name of vaccination scenarios.
#' @param rev_cov A logical variable that determines whether to revise the SIA
#' coverage data from Montagu and use precise age at vaccination.

#' @examples
#'   create_vaccine_coverage_routine_sia (
#'   vaccine_coverage_folder    = "vaccine_coverage/",
#'   vaccine_coverage_subfolder = "scenarios/",
#'   coverage_prefix            = "coverage",
#'   touchstone                 = "_201910gavi-5_",
#'   antigen                    = "measles-",
#'   scenario_name              = "campaign-only-bestcase",
#'   rev_cov                    = TRUE
#'   )
create_vaccine_coverage_routine_sia <- function (vaccine_coverage_folder    = "",
                                                 vaccine_coverage_subfolder = "",
                                                 coverage_prefix            = "",
                                                 touchstone                 = "",
                                                 antigen                    = "",
                                                 scenario_name              = "",
                                                 rev_cov
                                                 ) {

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
  vaccov <- fread (file = vaccine_coverage_file,
                   stringsAsFactors = F, na.strings = "<NA>")

  keep.cols <- c("vaccine", "country_code", "country", "year", "age_first", "age_last",
                 "age_range_verbatim", "target", "coverage")

  # select routine vaccination coverage
  routine <- vaccov [activity_type != "campaign", ..keep.cols]

  # only select campaigns with coverage > 0 and with information of target population size
  sia <- vaccov [activity_type == "campaign" & ( !is.na(target) & !is.na(coverage) & coverage != 0 ),
              ..keep.cols]

  # add columns for the number of reached and precise age at vaccination of SIAs
  sia [, `:=` (a0 = 0, a1 = 0,
               reached = as.numeric(target) * as.numeric(coverage))]

  if ((rev_cov) & (nrow(sia) > 0)) {

    # remove all whitespace from age_range_verbatim, put all in lowercase
    sia$age_range_verbatim <- tolower (gsub("[[:space:]]", "", sia$age_range_verbatim))

    # adjust format for following extraction of precise age
    sia [age_range_verbatim == "9-5y" , age_range_verbatim := "9m-5y"]
    sia [age_range_verbatim == "1to4", age_range_verbatim := "1-4y"]
    sia [age_range_verbatim == "<5y", age_range_verbatim := "0-<5y"]
    sia [age_range_verbatim == "<15y", age_range_verbatim := "0-<15y"]
    sia [age_range_verbatim == "6m+" , age_range_verbatim := "6m-100y"]
    sia [age_range_verbatim == "1yschool" , age_range_verbatim := "1-6y"]

    sia [, `:=`(age_first_verb = stringr::str_extract (age_range_verbatim, "[^-]+"),
                age_last_verb  = stringr::str_extract (age_range_verbatim, "[^-]+$"))]

    sia [, `:=`(age_first_num  = as.double (stringr::str_extract (age_first_verb, "\\d+")),
                age_first_unit = stringr::str_extract (age_first_verb, "\\D+$"),
                age_last_num   = as.double (stringr::str_extract (age_last_verb, "\\d+")),
                age_last_unit  = stringr::str_extract (age_last_verb, "\\D+$"))]

    sia [is.na(age_first_unit), age_first_unit := age_last_unit]
    sia [, `:=` (age_first = as.double (age_first), age_last = as.double (age_last))]

    # use precise age for those <3 years old
    sia [age_first_unit == "m" ,
         age_first := ifelse(age_first_num < 36, age_first_num/12, round(age_first_num/12))]
    sia [age_last_unit == "m" ,
         age_last := ifelse(age_last_num < 36, age_last_num/12, round(age_last_num/12))]

    # recalculate the size of target population
    get_targetpop <- function (icty, iyr, iage1, iage2){
      targetpop <- sum (data_pop [country_code == icty & year == as.integer(iyr) &
                                    age_from %in% (ceiling(iage1):floor(iage2)), value])
      if (ceiling(iage1) > iage1){
        targetpop <- targetpop +
          (ceiling(iage1)-iage1)*data_pop [country_code == icty & year == as.integer(iyr) &
                                     age_from == ceiling(iage1)-1, value] }

      if (iage2 > floor(iage2)){
        targetpop <- targetpop +
          (iage2-floor(iage2))*data_pop [country_code == icty & year == as.integer(iyr) &
                                           age_from == floor(iage2)+1, value]}

      return (targetpop)
    }

    sia [, target2 := get_targetpop(country_code, year, age_first, age_last),
         by = seq_len (nrow(sia))]

    # use VIMC results for special age ranges
    sia [age_range_verbatim %in% c("defaultageandgender", "school-age", "16-35yfemales", "17-24ymales"),
          target2 := target]

    sia [, coverage2 := reached/target2]
    sia [coverage2 > 1, coverage2 := 1]

    # remove temporary columns
    sia [, `:=` (target = target2, coverage = coverage2)]
    sia <- sia [, .SD, .SDcols = !c("age_first_verb", "age_last_verb",
                                    "age_first_num", "age_first_unit",
                                    "age_last_num", "age_last_unit",
                                    "target2", "coverage2")]
    }

  # -----------------------------------------------------------------------------------
  # calculate values for a0 and a1 to match the age groups
  # for age <= 3 years, age is in weeks
  # for age >= 3 years, (156 weeks for up to age 3) + (remaining years above 3)
  #      age 1 year ~ 52; 2 year ~ 104; 3 year ~ 156; 4 year ~ 157; 5 year ~ 158
  find.a <- function (x) {
    t0 <- ifelse (x <= 3, round (x * 52), round ((x - 3) + (3 * 52) ) )
    return (t0) }

  sia [, `:=` (a0 = find.a(age_first), a1 = find.a(age_last))]

  # set age == 0 year as 1 week
  sia [a0 == 0, a0 := 1]
  sia [a1 == 0, a1 := 1]

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
  # applies in both past and future years. *target* is the number of individuals
  # in the target population.
  #
  # This is always shown for campaigns, and is now specified
  # at a national level. (See ‘coverage’ section above.)
  #
  # For routine, target is shown as NA, which means you should assume the target
  # population matches the population shown in the demographic data downloads
  # for the corresponding ages (age_first and age_last).
  #
  # ----------------------------------------------------------------------------
  # REVISIONS FOR 2021 RUNS:
  # The coverage is revised by recalculating the number of reached by campaign
  # (target*coverage, numerator) and the size of targeted age group using the
  # UNWPP data (denominator). This is to improve the precision of coverage
  # estimates and age at vaccination.
  #
  # However, one limitation remains - for those with Montagu coverage = 1, the
  # number of reached equals to the size of targeted population. The original
  # number 'reached' cannot be obtained.
  # ----------------------------------------------------------------------------

  # write vaccine coverage data for routine vaccination and SIAs
  fwrite (x = routine, file = routine_coverage_file)
  fwrite (x = sia,     file = sia_coverage_file)

} # end of function -- create_vaccine_coverage_routine_sia
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Expand the matrix to a different dimension
#
#' This function returns an expanded contact matrix with the specified age
#' structure for model inputs.
#'
#' @param A A matrix to be expanded.
#' @param expand_rows A number of times to repeat each row.
#' @param expand_cols A number of times to repeat each column.
#' @param rescale_rows A logical variable to control whether to rescale the
#' expanded rows.
#' @param rescale_cols A logical variable to control whether to rescale the
#' expanded columns.
#' @return An expanded matrix with rescaling if applicable.
#' @examples
#' expandMatrix (matrix(1:9,3,3), 2, 1, FALSE, FALSE)
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
#' Create a vaccine coverage file for no-vaccination scenario
#'
#' A csv file for 0\% coverage vaccination scenario is created based on the csv
#' file of any of other vaccination scenarios. This is done because the coverage
#'  file downloaded from VIMC for the no-vaccination scenario is empty.
#'
#' @param no_vaccination_coverage_file File name given to the 0\% vaccine
#' coverage scenarios
#' @param vaccination_coverage_file File used to generate vaccine coverage for
#' 0\% scenarios
#' @examples
#' create_no_vaccination_coverage_file (
#'  no_vaccination_coverage_file = "no-vaccination.csv",
#'  vaccination_coverage_file    = "mcv1-default.csv"
#'  )
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
#' Create a vaccine coverage file for campaign-only scenario
#'
#' A csv file for campaign-only vaccination scenario is created based on the
#' coverage file for (routine + campaign) scenarios, by setting routine coverage
#'  to zero. This is done because routine coverage values are needed even if
#'  they are only zeros to run campaign-only vaccination scenarios.
#'
#' @param campaign_only_vaccination_coverage_file File name given to the
#'  campaign-only scenarios
#' @param routine_campaign_vaccination_coverage_file File used to generate
#'  vaccine coverage for campaign-only scenarios
#' @examples
#' create_campaign_vaccination_coverage_file (
#'  campaign_only_vaccination_coverage_file    = "campaign-only-bestcase.csv",
#'  routine_campaign_vaccination_coverage_file = "campaign-bestcase.csv"
#'  )
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
#' Tailor the data structure for life expectancy by age and year
#'
#' Tailor the \code{\link{data_lexp_remain}} data to the format for processing burden
#' estimates. Linear interpolation between calender years was applied.
#'
#' @param sel_countries ISO-3 codes for countries included for evaluation. If
#' "all", all countries in the original data are selected.
#' @examples
#' lexp_remain <- tailor_data_lexp_remain (sel_countries = c("AGO","BGD"))
tailor_data_lexp_remain <- function (sel_countries = "all"){

  lexp_remain <- setDT (data_lexp_remain)
  lexp_remain <- lexp_remain [year >= 1980]

  if (sel_countries[[1]] != "all") {
    lexp_remain <- lexp_remain [country_code %in% sel_countries]
  }

  # copy values for year 2095 to year 2100
  lexp_remain <- rbind (lexp_remain, copy (lexp_remain [year == 2095])[, year := 2100])

  # calculate the difference between years
  setorder (lexp_remain, country_code, age_from, year)
  lexp_remain [ , diffy := value - shift(value) , by = .(country_code, age_from)]

  # interpolate a linear trend for between-years
  lexp_remain_yr <- copy (lexp_remain) [year == 1980]
  for (btwyr in 0:4) {
    dt <- copy (lexp_remain) [year != 1980]
    dt [, `:=` (year = year - btwyr,
                value = value - btwyr*(diffy/5))]

    lexp_remain_yr <- rbind (lexp_remain_yr, dt)
  }

  # interpolate a linear trend for between-ages
  setorder (lexp_remain_yr, country_code, year, age_from)
  lexp_remain_yr [ , age_mean := (age_from + age_to)/2]

  lexp_remain_full <- lexp_remain_yr [, .(value = approx(age_mean, value, xout = 0:100)$y),
                                      by = .(country_code, country, year)]
  lexp_remain_full [, age := rep(0:100, length(unique(year))*length(unique(country_code)))]

  return (lexp_remain_full)

} # end of function -- tailor_data_lexp_remain
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Create variable inputs for probabilistic sensitivity analysis
#'
#' Generate latin hyper cube samples of input parameters for probabilistic
#' sensitivity analysis (PSA). Four variables are included in the output: (1)
#' "take1_input"-intercept of the regression model for age-varying vaccine
#' protection under 'take' assumption, (2) "take2_input"-coefficient of the
#' regression model for age-varying vaccine protection under 'take' assumption,
#' (3) "take3_input"-vaccine protection for the second dose under 'take'
#' assumption, (4) "mortality_input"-proportion of change in case fatality
#' ratios.
#'
#' @param psa An integer indicating the total runs of PSA.
#' @param seed_state An integer for selecting the seed number.
#' @param psadat_filename A file name containing ID and input parameters of each
#'  run of PSA.
#' @examples
#' create_PSA_data (psa = 200, seed_state = 1)
create_PSA_data <- function (psa             = 0,
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
  cube <- lhs::randomLHS (n = psa,
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

  psadat [, vaceffbyage_b := truncnorm::qtruncnorm (cube [, 1],
                                                    a    = 0,
                                                    # a    = (mean_slope - 3 * sd_slope), == -0.0004079801
                                                    b    = (mean_slope + 3 * sd_slope),
                                                    mean = mean_slope,
                                                    sd   = sd_slope
                                                    )
          ]

  # ----------------------------------------------------------------------------

  # vaccine efficacy (dose 2)
  # mean 98% -- +- ~ 2%
  psadat [, take3_input := truncnorm::qtruncnorm (cube [, 1],
                                                  a    = (0.98 - 0.02),
                                                  b    = (0.98 + 0.02),
                                                  mean = 0.98,
                                                  sd   = (0.02/3)
                                                  )
          ]


  # proportional change in case fatality rate
  # truncated lognormal distribution for up to 25% change
  psadat [, mortality_input := truncnorm::qtruncnorm (cube [, 2],
                                                      a    = (1 - 0.25),
                                                      b    = (1 + 0.25),
                                                      mean = 1,
                                                      sd   = (0.25/3)
                                                      )
          ]

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
    fwrite (x = psadat, file = psadat_filename)
  }

  # return psa data table
  return (psadat)

} # end of function -- create_PSA_data

# Codes for generating PSA variables in the model without age-dependent vaccine
# efficacy in take. (remvoed fron "runScenario")
# # create csv if file does not exist
# psa_var <- data.table(
#   run_id = c(1:psa),
#   take1_input     = rep (NA, psa),
#   take2_input     = rep (NA, psa),
#   take3_input     = rep (NA, psa),
#   mortality_input = rep (NA, psa)
# )
#
# for(t in 1:3) {
#
#   # use 5% difference in take
#   tk <- runif (n   = psa,
#                min = (take[t] * 0.95),
#                max = (take[t] * 1.05))
#
#   # vaccine efficacy is bounded between 0 and 1
#   tk [which(tk < 0)] <- 0
#   tk [which(tk > 1)] <- 1
#
#   psa_var [, paste0 ("take", t, "_input")] <- tk
# }
#
# # use 25% difference in mortality, will be multiplied by CFR in each country
# psa_var [, "mortality_input"] <- runif (n  = psa,
#                                         min = 0.75,
#                                         max = 1.25)
#
# fwrite (psa_var, paste0 ("input/", data_psa))
# ------------------------------------------------------------------------------


# Not in use the latest version, where deaths are estimated by "function/estimate_deaths_dalys.r"
# # ------------------------------------------------------------------------------
# #' Tailor the data structure for case fatality ratios (CFRs)
# #'
# #' Tailor \code{\link{data_cfr}} to the format for processing burden estimates.
# #'  Changes include calender year and age structures.
# #'
# #' @param sel_countries ISO-3 codes for countries included for evaluation. If
# #' "all", all countries in the original data are selected.
# #' @examples
# #' cfr.year.all <- tailor_data_cfr (sel_countries = "all")
# tailor_data_cfr <- function (sel_countries){
#
#   cfr <- setDT(data_cfr)
#
#   if (sel_countries[[1]] != "all") {
#     cfr <- cfr [Code %in% sel_countries]
#   }
#
#   # set cfrs before 2000 to 2000
#   cfr.year <- rbindlist(lapply(1980:1999, function(i) copy(cfr[Year ==2000,])[,Year := i]))
#   cfr <- rbind(cfr.year, cfr)
#
#   # set Kosovo mortaity to Serbia as it was once Serbia
#   cfr.xk <- rbindlist(lapply("XK", function(i) copy(cfr[Code == "SRB",])[,Code := i]))
#   cfr.xk[, Country := "Kosovo"]
#   cfr <- rbind(cfr, cfr.xk)
#
#   # expand cfr by age so that it uses the right value for <5 and 5-9 and 0 for >=10
#   cfr.year.all <- rbindlist(lapply(0:100, function(i) copy(cfr)[, age:=i]))
#   over10 <- T
#
#   if (over10){
#     cfr.year.all[, cfr.value :=  if (age < 5) under5 else if (age > 4 & age <10) over5 else 0,
#                  by = c("Code", "Year", "age") ]
#   } else {
#     cfr.year.all[, cfr.value :=  if (age < 5) under5 else if (age > 4 ) over5,
#                  by = c("Code", "Year", "age") ]
#   }
#   cfr.year.all$Year = as.integer(cfr.year.all$Year)
#
#   return(cfr.year.all)
# }
# # ------------------------------------------------------------------------------


# Not in use the latest version
# # ------------------------------------------------------------------------------
#  timelinessCov
#
#  get vaccination coverages according to the age groups used for model inputs.
#  @param tim_shape
#  @param target_age
#  @param max_cov
#  @examples
#  timelinessCov()
# timelinessCov <- function (tim_shape,
#                            target_age,
#                            max_cov) {
#
#   coverage_byage_cumulative <- c(
#     rep(0,(target_age-1)),
#     rep(max_cov,(254-target_age+1))
#   )
#
#   beta_function <- cumsum(
#     dbeta(
#       (1:52/52),
#       4,
#       tim_shape
#     )
#   )/max(
#     cumsum(
#       dbeta(
#         (1:52/52),
#         4,
#         tim_shape
#       )
#     )
#   )
#
#   coverage_byage_cumulative[39:(39+51)] <- max_cov*beta_function
#   coverage_byage=rep (0, 254)
#   coverage_byage [2:length(coverage_byage_cumulative)] <- (
#     coverage_byage_cumulative [2:length(coverage_byage_cumulative)]
#     -coverage_byage_cumulative [1:(length(coverage_byage_cumulative)-1)]
#   )/(
#     1 - coverage_byage_cumulative [1:(length(coverage_byage_cumulative)-1)]
#   )
#
#   coverage_byage [which (is.na (coverage_byage))] <- 1
#
#   return (coverage_byage)
#
# } # end of function -- timelinessCov
# # ------------------------------------------------------------------------------
