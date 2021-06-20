# dxplot.R
# Function for diagnosing model problems.


# ------------------------------------------------------------------------------
#' diagnostic_plots
#'
#' Create diagnostic plots of vaccine coverages and burden estimates (cases,
#' deaths, dalys) for selected countries.
#' @param vaccine_coverage_folder Folder name for the vaccine coverage files.
#' Include a slash in the end.
#' @param coverage_prefix Prefix of vaccine coverage file names.
#' @param touchstone Version note used by VIMC. Include a underscore at the
#' beginning and the end.
#' @param antigen Disease name used by VIMC.
#' @param scenarios Name of vaccination strategies included for evaluation.
#' @param base_scenario Name of the vaccination strategy that is used as the
#' baseline for calculating averted burden.
#' @param burden_estimate_folder A folder name for the file which contains the
#' model outputs for evaluation.Include a slash at the end.
#' @param plot_folder A folder name for the diagnostic plots to be saved.
#' @param group_name A modelling group name used by VIMC.
#' @param countries A vector of ISO country codes for plotting.
#' @param cfr_option The method used for specifying the trend of case fatality
#' ratios (CFRs): "Wolfson" or/and "Portnoy".
#' @param psa Number of runs for probabilistic sensitivity analysis.
#' @param start_year Beginning year of the plot.
#' @param end_year End year of the plot.
#' @param compare_plots A logical variable to determine whether to combine model
#'  outputs across different scenarios and compare them in the same plot.
#' @examples
#' diagnostic_plots (
#'   vaccine_coverage_folder    = "vaccine_coverage/",
#'   coverage_prefix            = "coverage",
#'   touchstone                 = "_201910gavi-5_",
#'   antigen                    = var$antigen,
#'   scenarios                  = scenarios[1:3],
#'   base_scenario              = base_scenario,
#'   burden_estimate_folder     = "central_burden_estimate/",
#'   plot_folder                = "plots/",
#'   group_name                 = "LSHTM-Jit-",
#'   countries                  = c("BGD","ETH"),
#'   cfr_options                = c("Wolfson", "Portnoy"),
#'   psa                        = 200,
#'   start_year                 = 2000,
#'   end_year                   = 2100,
#'   compare_plots              = FALSE)
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
        stringr::str_replace (string      = burden_estimate_file,
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
      if (countries[[1]] == "all") {    # (Han: consider using countries[1]=='all', so multiple countries can be selected)
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
          scale_x_continuous (breaks = scales::pretty_breaks ()) +
          geom_point () +
          labs (title = countrycode::countrycode (sourcevar   = country_iso3_code,
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

          p <- ggplot2::ggplot(burden_estimate [country == country_iso3_code],
                               aes(x = year, y = get(toplot))) +
            scale_x_continuous (breaks = scales::pretty_breaks ()) +
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) +
            labs (title = countrycode::countrycode (sourcevar   = country_iso3_code,
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
        plots <- ggpubr::ggarrange (plotlist = plot_list,
                                    ncol = 2, nrow = 2)

        # print plots
        print (ggpubr::annotate_figure (plots,
                                        top = ggpubr::text_grob (scenario_name,
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

          p <- ggplot2::ggplot(all_burden [country == country_iso3_code],
                               aes(x     = year,
                                   y     = get (toplot),
                                   group = scenario,
                                   color = factor (scenario))) +
            scale_x_continuous (breaks = scales::pretty_breaks ()) +
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) +
            labs (title = countrycode::countrycode (sourcevar = country_iso3_code,
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

          p <- ggplot2::ggplot(all_burden_compare [country == country_iso3_code],
                               aes(x     = year,
                                   y     = get (toplot),
                                   group = scenario,
                                   color = factor (scenario))) +
            scale_x_continuous (breaks = scales::pretty_breaks ()) +
            stat_summary (fun = sum, geom = "line") +
            ylab (toplot) +
            labs (title = countrycode::countrycode (sourcevar   = country_iso3_code,
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

          p <- ggplot2::ggplot(all_burden_compare [country == country_iso3_code],
                               aes(x     = scenario,
                                   y     = get (toplot),
                                   fill = factor (scenario))) +
            geom_col() +
            ylab (toplot) +
            labs (title = countrycode::countrycode (sourcevar   = country_iso3_code,
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

          p <- ggplot2::ggplot(all_burden_compare [country == country_iso3_code],
                               aes(x     = scenario,
                                   y     = get (toplot),
                                   fill = factor (scenario))) +
            geom_col() +
            ylab (toplot) +
            labs (title = countrycode::countrycode (sourcevar   = country_iso3_code,
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
