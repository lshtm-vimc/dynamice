# data.R
# Documentation of datasets


#' UNWPP population estimates
#'
#' A dataset containing the UNWPP population estimates -- World Population Prospects 2019.
#'
#' @format A data table with 27137760 observations of 8 variables.
#' \describe{
#'   \item{country_code_numeric}{Country code numeric}
#'   \item{country_code}{ISO3 country code}
#'   \item{country}{Country name}
#'   \item{age_from}{age from (start-age)}
#'   \item{age_to}{age to (end-age)}
#'   \item{year}{Year}
#'   \item{gender}{Gender}
#'   \item{value}{Population size}
#' }
#' @source \url{https://population.un.org/wpp/Download/Standard/Population/}
#' @source {VIMC}
"data_pop"


#' Contact matrix by age
#'
#' A dataset containing the contact matrix for age between 0 and 100,
#'  adapted from UK population in the POLYMOD study.
#'
#' @format A data table with 10302 observations of 102 variables.
#' \describe{
#'   \item{contact.age}{age between 0 and 100}
#'   \item{V1-V101}{age between 0 and 100}
#' }
#' @source {Mossong J, Hens N, Jit M, et al.
#' Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases.
#' PLoS Med 2008; 5(3): e74.}
"data_contact"


#' R0 by country
#'
#' A dataset containing the R0 of measles in 112 countries.
#'
#' @format A data table with 336 observations of 3 variables.
#' \describe{
#'   \item{contact.age}{age between 0 and 100}
#'   \item{V1:V101}{age between 0 and 100}
#' }
#' @source {Estimates provided by and Ferrari M and team}
"data_r0"


#' Timeliness data
#'
#' A dataset containing the timeliness data estimates and final coverage.
#'
#' @format A data table with 34496 observations of 4 variables.
#' \describe{
#'   \item{country_code}{ISO3 country code}
#'   \item{age}{age in weeks. \code{NA} represents the overall population. For
#'   0-1 year old (0-52 weeks),  \code{timeliness} is presented each week; for
#'   1-3 year old (53-156 weeks), it is presented the first week of each month.}
#'   \item{timeliness}{timeliness coverage according to age}
#'   \item{prop_final_cov}{final coverage of the population. Exist only when
#'   age is \code{NA}}
#' }
#' @source {Clark, A and Sanderson, C. Timing of children's vaccinations in 45
#' low-income and middle-income countries: an analysis of survey data. Lancet
#' 2009; 373(9674): 1543-1549.}
#' @source {Check with Kevin/Mark for how the csv is processed}
"data_timeliness"


#' Life expectancy at birth
#'
#' A dataset containing estimates of life expectancy at birth.
#'
#' @format A data table with 137712 observations of 8 variables.
#' \describe{
#'   \item{country_code_numeric}{Numeric country code}
#'   \item{country_code}{ISO3 country code}
#'   \item{country}{Country name}
#'   \item{age_from}{age from (start-age)}
#'   \item{age_to}{age to (end-age)}
#'   \item{year}{Year}
#'   \item{gender}{Gender}
#'   \item{LE0}{Life expectancy at birth (\code{age_from=age_to=0})}
#' }
#' @source \url{https://population.un.org/wpp/Download/Standard/Population/}
#' @source {VIMC}
"data_lexp0"


#' Life expectancy for a specific age group (remain)
#'
#' A dataset containing estimates of life expectancy for a specific age group.
#'
#' @format A data table with 601920 observations of 8 variables.
#' \describe{
#'   \item{country_code_numeric}{Numeric country code}
#'   \item{country_code}{ISO3 country code}
#'   \item{country}{Country name}
#'   \item{age_from}{age from (start-age)}
#'   \item{age_to}{age to (end-age)}
#'   \item{year}{Year}
#'   \item{gender}{Gender}
#'   \item{value}{Life expectacy for the age group between age_from and age_to}
#' }
#' @source \url{https://population.un.org/wpp/Download/Standard/Population/}
#' @source {VIMC}
"data_lexp_remain"


#' CFRs
#'
#' A dataset containing the case fatality ratios (CFRs) by calender year and age
#'  groups of under and over 5 years.
#'
#' @format A data table with 78477 observations of 7 variables.
#' \describe{
#'   \item{Country}{Country name}
#'   \item{Code}{ISO3 country code}
#'   \item{GBD Region}{GBD regions}
#'   \item{Grouping}{GBD regions further grouped by under-5 mortality rate,
#'   for example, LAC>=50 and LAC<50}
#'   \item{Year}{Year}
#'   \item{over5}{CFRs for the age group over 5 years old}
#'   \item{uner5}{CFRs for the age group under 5 years old}
#' }
#' @source {check with Kevin/Kaja/Mark}
"data_cfr"


#' CFRs for DALY calculation with Portnoy estimates
#'
#' A dataset containing the case fatality ratios (CFRs) by scenario, calender
#' year and age groups of under and over 5 years. using Portnoy estimates.
#'
#' @format A data table with 156240 observations of 45 variables.
#' \describe{
#'   \item{country_name}{Country name}
#'   \item{country_code}{ISO3 country code}
#'   \item{GBD Region}{GBD regions}
#'   \item{U5MR}{Groups for under-5 mortaltiy rate, >=50 or LAC<50}
#'   \item{year}{Year}
#'   \item{cfr_over5_scenario_s4}{CFRs for the age group over 5 years under one
#'   of the 10 coverage scenarios. The \code{s4} refers to declining trend after
#'    2018.}
#'   \item{cfr_over5_scenario_s6}{CFRs for the age group over 5 years under one
#'   of the 10 coverage scenarios. The \code{s6} refers to static trend of the
#'   2018 level in the future.}
#'   \item{cfr_uner5_scenario_s4}{CFRs for the age group under 5 years old,
#'   under one of the 10 coverage scenarios. The \code{s4} refers
#'   to declining trend after 2018.}
#'   \item{cfr_uner5_scenario_s6}{CFRs for the age group under 5 years old,
#'   under one of the 10 coverage scenarios. The \code{s6} refers to static
#'   trend of the 2018 level in the future.}
#' }
#' @source {Portnoy A, Jit M, Ferrari M, Hanson M, Brenzel L, Verguet S.
#' Estimates of case-fatality ratios of measles in low-income and middle-income
#' countries: a systematic review and modelling analysis. Lancet Glob Health
#' 2019; 7: e472–81.}
"data_cfr_portnoy"


#' CFRs for DALY calculation with Wolfson estimates
#'
#' A dataset containing the country-specific case fatality ratios among children
#'  under 5 years old, using Wolfson estimates.
#'
#' @format A data table with 294 observations of 3 variables.
#' \describe{
#'   \item{country_code}{ISO3 country code}
#'   \item{Country}{Country name}
#'   \item{CFR}{Case fatality ratios}
#' }
#' @source {Wolfson LJ, Grais RF, Luquero FJ, Birmingham ME, Strebel PM.
#' Estimates of measles case fatality ratios: a comprehensive review of
#' community-based studies. Int J Epidemiol 2009; 38: 192–205.}
"data_cfr_wolfson"


#' template for outputs
#'
#' A template dataset containing countries, years, ages, and other outputs.
#'
#' @format A data table with 8908200 observations of 9 variables.
#' \describe{
#'   \item{disease}{Disease (="Measles")}
#'   \item{year}{Year, from 2000 to 2100}
#'   \item{age}{Age, from 0 to 100 years old}
#'   \item{country}{ISO3 country code}
#'   \item{country_name}{Country name}
#'   \item{cohort_size}{Population size}
#'   \item{deaths}{Disease-specific deaths}
#'   \item{cases}{Disease-specific cases}
#'   \item{dalys}{Disease-specific disability-adjusted life years}
#' }
#' @source {VIMC}
"data_template"
