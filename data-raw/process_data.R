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
  data_contact_polymod <- local({

    # use POLYMOD matrix with physical contacts in Great Britain
    # obtain from the original paper:
    # https://doi.org/10.1371/journal.pmed.0050074.st005 (Table S8.4b)

    # data for all contacts (physical and non-physical) also available through 'socialmixr' R package
    # https://cran.r-project.org/web/packages/socialmixr/vignettes/introduction.html
    # contact_ori <- socialmixr::contact_matrix (polymod,
    #                                            countries = "United Kingdom",
    #                                            age.limits = seq(0,80,5),
    #                                            symmetric = FALSE)$matrix

    # ensure the matrix to show contactors in rows and contactees in columns
    contact_ori <- fread (paste0 (wd_rawdata, "polymod_physical_uk.csv"))
    contact_ori <- unname (t(contact_ori[, age := NULL]))

    # prepare a contact matrix with yearly age groups between 0-100
    # according to single-age populations for each country at reference year (2020)
    # an earlier version based on the UK population at survey year (2005): data_contact_uk.Rdata
    # assume contacts of 70+ include those between 70-100
    load (file = "data/data_pop.rda")
    pop_2020     <- data_pop [year == 2020]
    country_list <- unique(pop_2020$country_code)   # 195 countries

    contact_polymod <- sapply (country_list, function(icty) {

      pop_ref <- pop_2020 [country_code == icty, value]

      # expand contactees (break down in to single-year age groups by proportion)
      contact_col <- matrix (0, ncol = 101, nrow = dim(contact_ori)[1])

      # 0-69 years old
      for (icol in 1:14){
        pop_prp <- pop_ref[5*(icol-1)+(1:5)] / sum (pop_ref[5*(icol-1)+(1:5)])
        contactees <- rep (contact_ori[, icol], 5)
        contact_col[, 5*(icol-1)+(1:5)] <- sweep (matrix (contactees, ncol = 5), 2, pop_prp, "*")
      }

      # 70-100 years old
      pop_prp <- pop_ref[71:101] / sum (pop_ref[71:101])
      contact_col[, 71:101] <- sweep (matrix (rep(contact_ori[, 15], 100-70+1) , ncol = 100-70+1),
                                    2, pop_prp, "*")

      # expand contactors (apply the average to different groups)
      contactor_0to69   <- matrix (rep (contact_col[1:14,], each = 5), nrow = 14*5)
      contactor_70to100 <- matrix (rep (contact_col[15,], each = 100-70+1), nrow = 100-70+1)
      contact_101 <- rbind (contactor_0to69, contactor_70to100)


      # # adjust contact reciprocity using UNWPP data at survey year (2005)
      # contact_rec <- matrix(0, 101, 101)
      # for (i in 1:101){
      #   for (j in 1:101) {
      #     contact_rec [i, j] <- (contact_101[i, j] * pop_uk[i] +
      #                                contact_101[j, i] * pop_uk[j])/(2*pop_uk[i])
      #   }
      # }

      return(contact_101)
    },
    simplify = FALSE, USE.NAMES = TRUE)

    return(contact_polymod)
  })


  ## Synthetic contact matrices
  data_contact_syn <- local({

    # use updated synthetic contact matrices 2020 (acknowledge to Kiesha Prem)
    # download overall contact without any particular assumptions - "contact_all.rdata"
    # original data also available as a csv file
    # urlfile = "https://raw.githubusercontent.com/kieshaprem/synthetic-contact-matrices/master/output/syntheticcontactmatrices2020/synthetic_contacts_2020.csv"
    # mydata <- readr::read_csv (url (urlfile))

    load (url ("https://github.com/kieshaprem/synthetic-contact-matrices/blob/master/output/syntheticcontactmatrices2020/overall/contact_all.rdata?raw=true"))
    load (file = "data/data_pop.rda")

    # use Ethiopia's synthetic matrix for Somalia (missing)
    contact_all[["SOM"]] <- contact_all[["ETH"]]

    contact_syn <- sapply (names(contact_all), function(icty) {

      contact_ori <- contact_all[[icty]]

      # prepare a contact matrix with yearly age groups between 0-100
      # according to single-age populations at the survey year (2005)
      # assume contacts of 75+ include those between 75-100
      pop_cty  <- data_pop [country_code == icty & year == 2020, value]

      # expand contactees (break down in to single-year age groups by proportion)
      contact_col <- matrix (0, ncol = 101, nrow = dim(contact_ori)[1])

      # 0-74 years old
      for (icol in 1:15){
        pop_prp <- pop_cty [5*(icol-1)+(1:5)] / sum(pop_cty[5*(icol-1)+(1:5)])
        contactees <- rep(contact_ori[, icol], 5)
        contact_col [, 5*(icol-1)+(1:5)] <- sweep (matrix (contactees, ncol = 5), 2, pop_prp, "*")
      }

      # 75-100 year old
      pop_prp <- pop_cty[76:101] / sum(pop_cty[76:101])
      contact_col[, 76:101] <- sweep (matrix (rep(contact_ori[, 16], 100-75+1) , ncol = 100-75+1),
                                      2, pop_prp, "*")   # 75-100 years old

      # expand contactors (apply the average to different groups)
      contactor_0to74   <- matrix (rep(contact_col[1:15,], each = 5), nrow = 15*5)
      contactor_75to100 <- matrix (rep(contact_col[16,], each = 100-75+1), nrow = 100-75+1)
      contact_101 <- rbind (contactor_0to74, contactor_75to100)

      # # adjust contact reciprocity using UNWPP data at reference year (2020)
      # contact_rec <- matrix(0, 101, 101)
      # for (i in 1:101){
      #   for (j in 1:101) {
      #     contact_rec [i, j] <- (contact_101[i, j] * pop_ref[i] +
      #                              contact_101[j, i] * pop_ref[j])/(2*pop_ref[i])
      #   }
      # }

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


  usethis::use_data(data_contact_polymod, # data_contact_uk
                    data_contact_syn,
                    data_lexp0,
                    data_cfr_wolfson,
                    overwrite = TRUE)
  #### ---------------------------------------------------------------------------
}
