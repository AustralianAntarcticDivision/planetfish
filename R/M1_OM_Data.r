## Specify parameters in the Operating Model


##* the next step for this function is to separate it into

#' Get operating model data
#'
#' Construct a list containing the operating model specifications
#' @inheritParams om
#' @export
get_om_data <- function(){
	## Define empty lists for OM parameters
	ctrl <- list() 	# General control parameters of simulations
	om <- list() 		# Parameters for the operating model
	sampling 	<- list()			#  Parameters to control population sampling
	### Control Parameters of Simulations
	ctrl$pin_casal_assess <- 1		# Run function 'run_complete_casal_assessment': Run CASAL assessment (1=yes, 0=no)
	ctrl$pin_update_casal_data <- 0		# In function 'run_complete_casal_assessment': Update casal data (1=yes, 0=no), used for scenario testing
	ctrl$pin_TAC_finder <- 1		# In function 'run_complete_casal_assessment': Run TAC finder after CASAL assessment (1=yes, 0=no)
	ctrl$pin_TAC_finder_om <- 0		# Run om_projection to find real TAC in final assessment year (1=yes, 0=no)
	ctrl$Assyr <- seq(2010,2010,2)		# Years with assessments e.g. seq(2010,2044,2)
	# casal_path <- rootdir	# Full path (e.g. rootdir and subdirectory) for CASAL files (can be in separate folder)
	ctrl$casal_path <- ""
	ctrl$inputprefix <- "casal_"
	ctrl$outputprefix <- "casalout_"
	ctrl <- update_casal_file_names(ctrl)		# Name casal input files
	# TAC Decision rules to split future TAC: 0: Use catch ratios from assessment year, 1: Catch ratios specified in 'specific_ratios'
	ctrl$use_specific_ratios <- 1
	ctrl$ssb_level_deplet <- 0.2		# Reference SSB level for depletion
	ctrl$ssb_level_target <- 0.5		# Reference SSB level for target (e.g. 0.5 for toothfish, 0.75 for icefish)
	ctrl$ref_prob_deplet <- 0.1		# Refernce probability of SSB dropping below ssb_level_deplet
	ctrl$ref_prob_target <- "median"	# Reference measure of SSB that needs to be at/above ssb_level_target
														# Using the median is hardwired, ref_prob_target is not actually used
														# (only introduced for clarification)
	# TAC Finder control parameters
	ctrl$pin_TAC_finder_new_DR <- 0		# Apply alternative decision rule and choose the overall lowest TAC (in om projections)
	ctrl$TAC_Finder_tol <- 0.01
	ctrl$TAC_Finder_TAC_multiplier <- 1.1
	ctrl$TAC_Finder_iterations_main <- 20
	ctrl$TAC_Finder_set_seed <- 0
	ctrl$TAC_Finder_rand_seed <- 0
	# Used for system calls (e.g. to CASAL): 0 = Windows version, 1 = Linux/Unix
	ctrl$linux <- 0
	### OM Model parameters (min, max values)
	## Ages
	om$age <- c(1,30)			# Age (min,max)
	om$ages <- seq(om$age[1], om$age[2])
	om$n_ages <- length(om$ages)
  ## years
	om$year <- c(1990,2010)		# Year (first,last). Make first year with 0 catch
	om$years <- seq(om$year[1], om$year[2])
	om$n_years <- length(om$years)
  ## Seasons
	om$season <- c(1,1)			# Seasons (min,max). Change also natM_props and growth_props
	om$seasons <- seq(om$season[1], om$season[2])
	om$n_seasons <- length(om$seasons)
  # Proportion of natM in each season
	om$natM_props <- c(1)
	# Prop of annual growth that has occurred by start each time step for calculation of length & weights etc
	om$growth_props <- c(0.5)
	## Regions
	om$region <- c(1,3)
	om$regions <- seq(om$region[1], om$region[2])
	om$n_regions <- length(om$regions)
	# Splitting recruitment to areas
	om$rec_area <- c(1,0,0)
	# Females first, Change: rec_sex, growth, WL, maturity, effsplitunit
	# Units or categories; If not sex, then change definition of ssb0 and ssb.current, and assessment input files (where with sex partition)!
	om$sex <- c("f","m")
	om$rec_sex <- c(0.5, 0.5) 		# Splitting recruitment to sex (was recunit)
	om$eff_sex <- c(0.5, 0.5) 		# Splitting fishing effort between sex (usually 0.5 & 0.5, ie F & M are targeted equally)
	### Fish population parametes
	# Stock recruitment relationshipRecruitment
	om$pin_rec <- "BH"
	# recruitment mu and sigma
	##* these are specified in an unusual manner
	om$rec_mu <- 5000000
	om$rec_sigma <- 0.6
	om$rec_h <- 0.7
	## this is an input into the SR function
	om$rec_series <- 0 ##* does this specify whether historic recruitment is used?
	# Growth by Von Bertalanffy function (Linf, K, t0, CV)
	om$growth <- list()
	om$growth[[om$sex[1]]] <- c(2870.8, 0.02056, -4.2897, 0.1)
	om$growth[[om$sex[2]]] <- c(2870.8, 0.02056, -4.2897, 0.1)
	# Weight-Length relationship (WLa, WLb)  (mm to tonnes)
	om$WL <- list()
	om$WL[[om$sex[1]]] <- list(a=0.00000000000259, b=3.2064)
	om$WL[[om$sex[2]]] <- list(a=0.00000000000259, b=3.2064)
	# Maturity: "logistic" (x50, x95), "ramp" (start, peak), provide" (provide a vector with all values)
	om$pin_mat <- "ramp"
	om$maturity <- list()
	om$maturity[[om$sex[1]]] <- list(start=11, peak=17)
	om$maturity[[om$sex[2]]] <- list(start=11, peak=17)
	# Natural mortality
	om$natM <- rep(0.155, om$n_ages)		# for testing: om$natM  <- seq(0.005,0.15,0.005)
	# initial age structure
	om$pin_init_age <- "stoch"
	#om$init_age_series <- rep(1e8, om$n_ages) # uses for pin_init_age = "provide" and "user"
	# Choose method to calculate SSB0: "stoch" (stochastic), "init_age" (???), "casal" (as in CASAL)
	om$B0_calculation_method <- "casal"
	om$niter_ssb <- 101	# If pin = "stoch": Number of iterations
	## Movement
	om$move_type <- "dframe"
	om$names_ages <- paste0("age_", om$ages)
	## move rates by age
	##* note this isn't retained in 'om'
	move_by_age <- as.data.frame(matrix(round(c(ogive("logistic", om$ages, c("x50"=5, "x95"=5))*0.7,
	                                         ogive("logistic", om$ages, c("x50"=10, "x95"=10))*0.7),4),
	                                 nrow=2, ncol=om$n_ages, byrow=TRUE))
	## dataframe of the movement rules
	om$move_rules <- data.frame("Origin" = c(1,2),"Destination" = c(2,3), "Sex" = c(0,0),
  	                          "Year" = c(0,0), "Season" = c(0,0), setNames(move_by_age, om$names_ages))
	## separate dataframe for the movement of tags
	om$move_tag_rules <- om$move_rules
	### Fishery parameters
	om$fishery <- c("Trawl","LL1","LL2")	#, "IUU")		# Fisheries or gear types
	om$n_fisheries <- length(om$fishery)
	# Selectivity specification pin
	om$pin_sel <- list()
	om$pin_sel[[om$fishery[1]]] <- "dbnormal"
	om$pin_sel[[om$fishery[2]]] <- "dbnormal"
	om$pin_sel[[om$fishery[3]]] <- "dbnormal"
	## selectivity parameters
	om$select <- list()
	# params[1]: a = top
	# params[2]: sigma.left
	# params[3]: sigma.right
	om$select[[om$fishery[1]]] <- list(top=5, sigma_left=2, sigma_right=5)
	om$select[[om$fishery[2]]] <- list(top=10, sigma_left=3, sigma_right=8)
	om$select[[om$fishery[3]]] <- list(top=15, sigma_left=6, sigma_right=8)
	## max harvest rate
	om$max_harvest <- 0.95
	### Initial / input catch and effort
	# Initial 'catch' and 'effort', for all years, fishery, seasons and regions
	# om$catch (intended catch) is stored in mod$real_catch (actual catch taken)
	# 	which may be lower than om$catch if the catch cannot be taken
	# When a TAC is available (eg from an assessment), TAC value will be used
	# 	(after applying an implementation error) and written into om$catch
	om$catch <- create_array(om, sampling, type="fishery")
	## define the catches for each fishery
	om$catch[,"Trawl",1,1] <- rep(1000, om$n_years)
	om$catch[,"LL1",1,2] <- rep(1000, om$n_years)
	om$catch[,"LL2",1,3] <- rep(1000, om$n_years)
	om$catch[1,,,] <- 0		# Catch in first year set to 0, such that SSB in first year is unfished biomass
	# Effort: values are currently not used
	om$effort <- om$catch
	# When TAC is calculated in assessment: Define here the actual catch splits between fisheries
	om$catchsplits <- create_array(om, sampling, type="fishery")
	## Define the catch splits for each fishery
	om$catchsplits[,"Trawl", 1, 1] <- rep(0.2, om$n_years)	# Trawl		# Should sum up to 1 within a year (across fishery, season & region)
	om$catchsplits[,"LL1", 1, 2] <- rep(0.4, om$n_years)	# LL1
	om$catchsplits[,"LL2", 1, 3] <- rep(0.4, om$n_years)	# LL2
	# Implementation and Observation(reporting) errors for catch
	om$catch_obs_error <- 0		# SD of lognormal distribution for difference between real catch and observed catch
	om$TAC_implement_error <- 0		# SD of lognormal distribution for difference between TAC and real catch taken
	## OM Projections
	om$y_proj <- 35		# Years that OM is projected for decision rules (similarly to CASAL projections)
	om$y_proj_alt <- 5		# Years that OM is projected for alternative decision rules
	om$N_proj <- 500		# Number of projections
	### Sampling Parameters: All parameters are by area, year and season
	## Is there sampling, either TRUE: Without sampling or FALSE: With sampling
	sampling$pin_sampling <- TRUE
	## Length classes for length samples
	sampling$len_classes <- seq(100,2000,50)
	sampling$n_lengths <- length(sampling$len_classes)
	## Catch at age/length: These numbers are for each year, season and region (if catch > 0)
	sampling$catchage_N <- 1000			# Total sample size for catch age (equally divided into sex, becomes catchsample.N in fleet, may need to be catch-weighted in sex-aggregated model!)
	sampling$catchlen_N <- 1000			# Total sample size for catch length (equally divided into sex, becomes catchsample.N in fleetlen, may need to be catch-weighted in sex-aggregated model!)
	# Ageing error - yet to be implemented
	# need to add ageing error to 'sample_survey'
	sampling$ageing_error <- 0.0
	## Survey index of absolute numbers at age/length
	sampling$survey_q <- 0.1			# Survey catchability
	sampling$survey_var <- 0.0			# Variance to add error to overall survey index
	sampling$survey_var_N <- 0.0			# Variance to add error to overall numbers (the same variance for all ages/lengths)
	sampling$survey_var_n <- 0.0			# Variance to add error to individual ages/lengths
	sampling$survey_var_units <- "lognormal"	# Error distribution of abundance index in units (normal, lognormal) (used to create abundance index)
	## CPUE as Index of relative abundance
	sampling$cpue_q <- 0.0001 		# Catchability for cpue ##* this needs to be by fishery
	sampling$cpue_var_real <- 0.2 			# Variance to create abundance index
	sampling$cpue_var_real_units <- "lognormal"	# Error distribution of abundance index in units (normal, lognormal) (used to create abundance index)
	### Tagging: total numbers for both sexes
	## Is tagging used either FALSE (Without tagging) or TRUE (With tagging)
	sampling$pin_tagging <- TRUE
	## sampling method either Pooled ("Pool") or Individual ("Ind")
	sampling$pin_tag_Method <- "Pool"
	## Random tagging FALSE: Without random tagging or TRUE: With random tagging
	sampling$pin_tag_Random <- TRUE
	# Specify how tags are released; "Fixed" Fixed number of tags, "Catch" Depending on catch
	sampling$pin_tag_N <- "Fixed"
	# For pin_tagging_N = 1: Total numbers of tags by a fishery, year, season and areas (if catch > 0)
	sampling$tag_N <- c(100, 100, 100)
	# For pin_tagging_N = 2: Total numbers of tags per tonne of fish caught (by sex)
	sampling$tag_rate <- c(2, 2, 2)
	## Tagging selectivity
	sampling$pin_tag_sel <- list()
	sampling$pin_tag_sel[[om$fishery[1]]] <- om$pin_sel[[om$fishery[1]]]
	sampling$pin_tag_sel[[om$fishery[2]]] <- om$pin_sel[[om$fishery[2]]]
	sampling$pin_tag_sel[[om$fishery[3]]] <- om$pin_sel[[om$fishery[3]]]
	## define the selectivity parameters
	sampling$tag_select <- list()
	sampling$tag_select[[om$fishery[1]]] <- om$select[[om$fishery[1]]] # c(10,2,10)
	sampling$tag_select[[om$fishery[2]]] <- om$select[[om$fishery[2]]] # c(5,2,6)
	sampling$tag_select[[om$fishery[3]]] <- om$select[[om$fishery[3]]] # c(5,2,6)
	# In OM: Proportion of tagged fish that are removed immediately after tagging
	sampling$tag_mort <- 0.1
	sampling$tag_shedding <- 0.0084	# In OM: Ongoing tag loss/shedding rate
	## create a list to store the three parameter lists
	para <- list("control" = ctrl,
	             "om" = om,
	             "sampling" = sampling)
	## return the parameter list
	para
}
