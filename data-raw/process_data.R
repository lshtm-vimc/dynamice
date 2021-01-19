# process_data.R
# This file processes raw data and crate exported data for the package.

Process_data <- function(){

  wd_rawdata <- paste0(getwd(),"/data-raw/")

  require("data.table")

  #### data without further processing ------------------------------------------------------
  # population, timeliness, contact matrix, life expectancy
  data_pop		  	  <- fread (paste0 (wd_rawdata, "201910gavi-5_dds-201910_2_int_pop_both.csv"))
  data_timeliness   <- fread (paste0 (wd_rawdata, "timeliness_new.csv"))
  data_cfr_portnoy  <- fread (paste0 (wd_rawdata, "cfr_scenarios.csv"))
  data_contact      <- fread (paste0 (wd_rawdata, "uk_polymod_physical_101.csv"))  # use newly rescaled polymod to yearly age bands to 100
  data_r0			    	<- fread (paste0 (wd_rawdata, "r0.csv"))
  data_cfr          <- fread (paste0 (wd_rawdata, "cfrs_new_noage.csv"))
  data_lexp_remain  <- fread (paste0 (wd_rawdata, "201910gavi-5_dds-201910_2_life_ex_both.csv"))
  data_template	  	<- fread (paste0 (wd_rawdata, "central-burden-template.201910gavi-5.Measles_LSHTM-Jit_standard.csv"))
  data_psa200   	  <- fread (paste0 (wd_rawdata, "psa_variables.csv"))

  usethis::use_data(data_pop,
                    data_timeliness,
                    data_cfr_portnoy,
                    data_contact,
                    data_r0,
                    data_cfr,
                    data_lexp_remain,
                    data_template,
                    data_psa200)

  #### --------------------------------------------------------------------------------------

  #### data need initial processing  --------------------------------------------------------
  # Life expectancy at birth
  data_lexp0 <- local({

    lexp0_in <- fread (paste0 (wd_rawdata, "201910gavi-5_dds-201910_lx0_both.csv"))

    # Generate data of 2100
    lexp0    <- rbindlist(lapply(2100, function(i) copy(lexp0_in[year == 2099])[, year:=i]))
    lexp0    <- rbind(lexp0_in, lexp0)

    setorder(lexp0, country_code, year)
    colnames(lexp0)[8] <- "LE0"

    return(lexp0)
   })


  # CFRs for DALYs calculation - Wolfson
  data_cfr_wolfson <- local({

    cfr		<- fread (paste0 (wd_rawdata, "cfr.csv"))

    # change CFR percentages to rates between 0 and 1
    cfr [, CFR := CFR/100]

    return(cfr)
  })

  usethis::use_data(data_lexp0,
                    data_lexp_remain,
                    data_cfr_wolfson,
                    overwrite = TRUE)
  #### ---------------------------------------------------------------------------
}
