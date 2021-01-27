# process_data.R
# This file processes raw data and crate exported data for the package.

Process_data <- function(){

  wd_rawdata <- paste0(getwd(),"/data-raw/")

  #### data without further processing -----------------------------------------
  data_pop		  	  <- fread (paste0 (wd_rawdata, "201910gavi-5_dds-201910_2_int_pop_both.csv"))
  data_timeliness   <- fread (paste0 (wd_rawdata, "timeliness_new.csv"))
  data_cfr_portnoy  <- fread (paste0 (wd_rawdata, "cfr_scenarios.csv"))
  data_r0			    	<- fread (paste0 (wd_rawdata, "r0.csv"))
  data_cfr          <- fread (paste0 (wd_rawdata, "cfrs_new_noage.csv"))
  data_lexp_remain  <- fread (paste0 (wd_rawdata, "201910gavi-5_dds-201910_2_life_ex_both.csv"))
  data_template	  	<- fread (paste0 (wd_rawdata, "central-burden-template.201910gavi-5.Measles_LSHTM-Jit_standard.csv"))
  data_psa200   	  <- fread (paste0 (wd_rawdata, "psa_variables.csv"))

  usethis::use_data(data_pop,
                    data_timeliness,
                    data_cfr_portnoy,
                    data_r0,
                    data_cfr,
                    data_lexp_remain,
                    data_template,
                    data_psa200)

  #### -------------------------------------------------------------------------


  #### data need initial processing  -------------------------------------------
  ## UK contact matrix
  data_contact_uk <- local({

    # use rescaled POLYMOD matrix for yearly age bands between 0-100
    contact    <- fread (paste0 (wd_rawdata, "uk_polymod_physical_101.csv"))
    contact    <- contact [ , -"contact.age"]

    return(contact)
  })

  ## Synthetic contact matrices
  data_contact_syn <- local({

    # use updated synthetic contact matrices 2020 (acknowledge to Kiesha Prem)
    # download overall contact without any particular assumptions
    # rescaled matrices for yearly age bands between 0-100
    load(url("https://github.com/kieshaprem/synthetic-contact-matrices/blob/master/output/syntheticcontactmatrices2020/overall/contact_all.rdata?raw=true"))

    contact_syn <- sapply (names(contact_all), function(cty){
      # rescale for a country-specific matrix
      contact_ori <- contact_all[[cty]]
      contact_101 <- matrix (0, ncol = 101, nrow = 101)

      nagegrp     <- dim(contact_ori)[1]
      matvals     <- numeric(0)
      for (icol in 1:nagegrp){
        matvals   <- c(matvals, rep (contact_ori[,icol], 5, each = 5))
      }
      contact_101 [1:(nagegrp*5), 1:(nagegrp*5)] <- matrix (matvals,
                                                            nrow = nagegrp*5,
                                                            ncol = nagegrp*5)
      # contacts> 80 y/o set to previously
      contact_101[81:101, 81:101] <- contact_101[80, 80]
      contact_101[81:101, 1:80]   <- matrix (rep (contact_101[80, 1:80], 101-80), nrow = 101-80, byrow = TRUE)
      contact_101[1:80, 81:101]   <- matrix (rep (contact_101[1:80, 80], 101-80), ncol = 101-80)

      return(contact_101)
      },
      simplify = FALSE, USE.NAMES = TRUE)

    return(contact_syn)
  })

  ## Life expectancy at birth
  data_lexp0 <- local({

    lexp0_in <- fread (paste0 (wd_rawdata, "201910gavi-5_dds-201910_lx0_both.csv"))

    # Generate data of 2100
    lexp0    <- rbindlist(lapply(2100, function(i) copy(lexp0_in[year == 2099])[, year:=i]))
    lexp0    <- rbind(lexp0_in, lexp0)

    setorder(lexp0, country_code, year)
    colnames(lexp0)[8] <- "LE0"

    return(lexp0)
   })


  ## CFRs for DALYs calculation - Wolfson
  data_cfr_wolfson <- local({

    cfr		<- fread (paste0 (wd_rawdata, "cfr.csv"))

    # change CFR percentages to rates between 0 and 1
    cfr [, CFR := CFR/100]

    return(cfr)
  })

  usethis::use_data(data_contact_uk,
                    data_contact_syn,
                    data_lexp0,
                    data_cfr_wolfson,
                    overwrite = TRUE)
  #### ---------------------------------------------------------------------------
}
