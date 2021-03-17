# process_data.R
# This file processes raw data and crate exported data for the package.

Process_data <- function() {

  wd_rawdata <- paste0 (getwd(),"/data-raw/")

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

    # use POLYMOD matrix with physical contacts for the Great Britain
    # original data from https://doi.org/10.1371/journal.pmed.0050074.st005 (Table S5b)
    # rescaled matrices for yearly age bands between 0-100 years old

    contact     <- fread (paste0 (wd_rawdata, "polymod_physical_uk.csv"))
    contact_ori <- t (as.matrix (contact [ , -"age"]))  # transpose to make participants/contactors presented in rows)
    contact_101 <- matrix (0, ncol = 101, nrow = 101)

    nagegrp     <- dim(contact_ori)[1]  # 15 groups
    expd_age1 <- 5         # 0-70 years old: expand from 5-year bands
    expd_age2 <- 101-70    # 70-100 years old: expand from the oldest group (70+)

    for (icol in 1:nagegrp){
      contactees <- c(rep (contact_ori[1:(nagegrp-1),icol]/expd_age1, each = expd_age1),
                      rep (contact_ori[nagegrp, icol]     /expd_age2, each = expd_age2))
      contact_101 [,(icol-1)*expd_age1+(1:expd_age1)] <-
        matrix (rep (contactees, expd_age1), ncol = expd_age1)
    }
    contact_101 [,(nagegrp*expd_age1+1):101] <-
      matrix (rep (contact_101 [, nagegrp*expd_age1], expd_age2-expd_age1), ncol = expd_age2-expd_age1)

    return(contact_101)
  })


  ## Synthetic contact matrices
  data_contact_syn <- local({

    # use updated synthetic contact matrices 2020 (acknowledge to Kiesha Prem)
    # download overall contact without any particular assumptions - "contact_all.rdata"
    # rescaled matrices for yearly age bands between 0-100 years old

    load (url ("https://github.com/kieshaprem/synthetic-contact-matrices/blob/master/output/syntheticcontactmatrices2020/overall/contact_all.rdata?raw=true"))

    contact_syn <- sapply (names(contact_all), function(icty) {

      contact_ori <- contact_all[[icty]]
      contact_101 <- matrix (0, ncol = 101, nrow = 101)

      nagegrp     <- dim(contact_ori)[1]  # 16 groups
      expd_age1 <- 5         # 0-74 years old: expand from 5-year bands
      expd_age2 <- 101-75    # 75-100 years old: expand from the oldest group (75+)

      for (icol in 1:nagegrp){
        contactees <- c(rep (contact_ori[1:(nagegrp-1),icol]/expd_age1, each = expd_age1),
                        rep (contact_ori[nagegrp, icol]     /expd_age2, each = expd_age2))
        contact_101 [,(icol-1)*expd_age1+(1:expd_age1)] <-
          matrix (rep (contactees, expd_age1), ncol = expd_age1)
      }
      contact_101 [,(nagegrp*expd_age1+1):101] <-
        matrix (rep (contact_101 [, nagegrp*expd_age1], expd_age2-expd_age1), ncol = expd_age2-expd_age1)

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
