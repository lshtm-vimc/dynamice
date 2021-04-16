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

    # use POLYMOD matrix with all contacts in Great Britain
    # obtain data from 'socialmixr' R package
    # ensure matrix to show contactors in rows and contactees in columns
    # have briefly check the consistency with original data from
    # https://doi.org/10.1371/journal.pmed.0050074.st005 (Table S8.4a)

    library('socialmixr')
    contact_ori <- socialmixr::contact_matrix (polymod,
                                               countries = "United Kingdom",
                                               age.limits = seq(0,80,5),
                                               symmetric = FALSE)$matrix

    # there are no survey participants over 80 years old
    # assume 80-84, 85-89, 90-94, 95-99 to have same contacts as 75-79
    contact_5y <- unname (contact_ori[1:16, 1:16])                                     # 16x16
    contact_5y <- cbind (contact_5y, matrix(rep(contact_5y[,16], 4), ncol = 4))        # 16x20
    contact_5y <- rbind (contact_5y, matrix(rep(contact_5y[16,], each = 4), nrow = 4)) # 20x20

    # expand matrices from 5-year to 1-year age bands
    contact_101 <- matrix (0, ncol = 101, nrow = 101)
    for (icol in 1:20){
      contactees <- rep(contact_5y[1:20, icol]/5, each = 5)    # 0-99 years old
      contact_101 [1:100, 5*(icol-1)+(1:5)] <- matrix (rep(contactees, 5), ncol = 5)
      rm(contactees)
    }
    contact_101 [, 101] <- contact_101 [, 100]   # 100 years old
    contact_101 [101, ] <- contact_101 [100, ]

    # adjust contact reciprocity using UNWPP data at survey year (2005)
    load (file = "data/data_pop.rda")
    pop_uk      <- data_pop [country == "United Kingdom" & year == 2005, value]
    contact_rec <- matrix(0, 101, 101)
    for (i in 1:101){
      for (j in 1:101) {
        contact_rec [i, j] <- (contact_101[i, j] * pop_uk[i] +
                                   contact_101[j, i] * pop_uk[j])/(2*pop_uk[i])
      }
    }

    return(contact_rec)
  })


  ## Synthetic contact matrices
  data_contact_syn <- local({

    # use updated synthetic contact matrices 2020 (acknowledge to Kiesha Prem)
    # download overall contact without any particular assumptions - "contact_all.rdata"
    # original data is also available as a csv file
    # urlfile = "https://raw.githubusercontent.com/kieshaprem/synthetic-contact-matrices/master/output/syntheticcontactmatrices2020/synthetic_contacts_2020.csv"
    # mydata <- readr::read_csv (url (urlfile))

    load (url ("https://github.com/kieshaprem/synthetic-contact-matrices/blob/master/output/syntheticcontactmatrices2020/overall/contact_all.rdata?raw=true"))
    load (file = "data/data_pop.rda")

    contact_syn <- sapply (names(contact_all), function(icty) {

      # there are no data for contacts over 80 years old
      # assume 80-84, 85-89, 90-94, 95-99 to have same contacts as 75-79
      contact_5y <- contact_all[[icty]]                                                  # 16x16
      contact_5y <- cbind (contact_5y, matrix(rep(contact_5y[,16], 4), ncol = 4))        # 16x20
      contact_5y <- rbind (contact_5y, matrix(rep(contact_5y[16,], each = 4), nrow = 4)) # 20x20

      # expand matrices from 5-year to 1-year age bands
      contact_101 <- matrix (0, ncol = 101, nrow = 101)
      for (icol in 1:20){
        contactees <- rep(contact_5y[1:20, icol]/5, each = 5)    # 0-99 years old
        contact_101 [1:100, 5*(icol-1)+(1:5)] <- matrix (rep(contactees, 5), ncol = 5)
        rm(contactees)
      }
      contact_101 [, 101] <- contact_101 [, 100]   # 100 years old
      contact_101 [101, ] <- contact_101 [100, ]

      # adjust contact reciprocity using UNWPP data at reference year (2020)
      pop_ref     <- data_pop [country_code == icty & year == 2020, value]
      contact_rec <- matrix(0, 101, 101)
      for (i in 1:101){
        for (j in 1:101) {
          contact_rec [i, j] <- (contact_101[i, j] * pop_ref[i] +
                                   contact_101[j, i] * pop_ref[j])/(2*pop_ref[i])
        }
      }

      return(contact_rec)
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
