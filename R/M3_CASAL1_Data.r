#### Input Parameters used for CASAL Assessment (population.csl, estimation.csl and output.csl)
#### Returns 'datass'

#' Get CASAL parameters
#'
#' Get CASAL parameters
#' @inheritParams om
#' @export
get_casal_para <- function(para) {
  # extract the om and sampling parameters
	om <- para[["om"]]
	sampling <- para[["sampling"]]
	### Create file with all input data (to be read into CASAL population file)
	datass <- list()
	#### Define general model parameters
	Yr_last <- om$year[2]		# Default assessment year is final year of om period (eg. 2010)
	datass$year <- c(om$year[1],Yr_last,which(om$years==Yr_last))	# first, current, which year in om
	datass$years <- seq(datass$year[1],datass$year[2])				# Sequence of all years

	yrs <- 1:datass$year[3]									# Used to define years in some arrays below
	datass$age <- om$age
	datass$ages <- om$ages
	datass$season <- c(1,1,1)		# om$season: First season for growth & ageing, second season for fishery
	datass$seasons <- c(1)		    # om$seasons (time steps)
	# Fisheries: Choose the match between the fisheries in the OM and in the assessment
	# Each fishery is unique by selectivity, season and region
	datass$Fish <- matrix(c("Trawl", 1, "R1", "SelTrawl", "qTrawl", 0.25, 100, 100,
	                        "LL",    1, "R2", "SelLL", "qLL", 0.75, 100, 100),
	                      ncol=8, byrow=TRUE, dimnames=list(c(),c("Fishery","Season","Region",
	                      "Sel","q","ProjCatch","catchage_N", "catchlen_N")))
	datass$match_fishery <- matrix(c("Trawl","LL1","LL2","Trawl","LL", "LL"),
	                               ncol=3, byrow=TRUE, dimnames=list(c("OM","Ass"),c()))
											# How are fisheries from the OM represented in the Assessment, e.g.
											# three OM fisheries ("Trawl","LL1","LL2"), but only two Ass fisheries ("Trawl","LL"):
											# matrix(c("Trawl","LL1","LL2",
											# 		   "Trawl","LL", "LL"), ncol=3, byrow=TRUE, dimnames=list(c("OM","Ass"),c()))
	# Check
	###*** replace with warning()
	if(sum(!om$fishery %in% datass$match_fishery["OM",])>0) print("Mismatch in OM fisheries (datass$match_fishery)")	# Warning
	# Regions: relationship between regions in OM and Assessment (i.e. allow for region mismatch between OM and assessment):
	datass$match_region <- matrix(c(  1,   2,   3, "R1","R2","R2"), ncol=3,
	                              byrow=TRUE, dimnames=list(c("OM","Ass"),c()))
											# How are regions from the OM represented in the Assessment, e.g.
											# three OM Regions (1,2,3), but only two Ass region (1,2):
											# matrix(c(1, 2, 3,
											#  		   1, 1, 2), ncol=3, byrow=TRUE, dimnames=list(c("OM","Ass"),c()))
	datass$regions <- sort(unique(datass$match_region["Ass",]))
	# Check
	###*** replace with warning()
	if(sum(!om$regions %in% datass$match_region["OM",])>0) print("Mismatch in OM regions (datass$match_region)")	# Warning
	datass$n_stocks <- 1			# Do not change
	# datass$stock_names			<- "Stock1"
	datass$y_proj <- 35			# Years for projection
	datass$by_sex	 <- 0			# 0: no sex partition, 1: with F & M
													# Sex partition not implemented
	datass$age_plus_group <- "True"				# Plus group for age data
	datass$len_plus_group <- "False"				# Plus group for length data
	datass$size_based <- "False"
	datass$class_mins <- sampling$len_classes # Must use len_classes (only purpose of len_classes in OM is to create samples)
	if (datass$len_plus_group == "False") 					# Additional value defines the upper limit of the last size class
			datass$class_mins <- c(datass$class_mins,(max(datass$class_mins)+(datass$class_mins[2]-datass$class_mins[1])))
	datass$sex_partition <- if (datass$by_sex==0) "False" else "True"
	datass$mature_partition <- "False"
	datass$sex <- c("f","m")			# Used to setup some FL objects within datass that are sex-based (observations)
	datass$initialB0 <- c(50000, 100000, 50000, 250000)		# B0: min/max bounds for starting values, min/max bounds for estimation
	datass$B0 <- runif(1, datass$initialB0[1],datass$initialB0[2])		# Initialisation
	# datass$B0						<- as.vector(round(mod$ssb0_init,0))		# If B0 is not estimamted
	###*** think about moving these comments somewhere into the doc files
	## CASAL sequence within a time step:
	## 		Ageing (in an age-based model)
	## 		Recruitment
	## 		Maturation (if maturity is a character in the partition)
	## 		Migration (if the model includes more than one area).
	## 		Growth (in a size-based model).
	## 		Natural and fishing mortality.
	## 		Disease mortality.
	## 		Tag release events.
	## 		Tag shedding rate.
	## 		Semelparous mortality.
	# Time steps in a year
	datass$time_steps <- datass$season[3]
	# Recruitment area: Recruitment in Casal can occur only to one area per stock
	#	Should recruitment occur to more areas: Recruitment into first area in list, then movement of proportion of recruits to all other areas
	datass$recruitment_areas <- "R1"					# Recruitment occurs in this area (only 1 area per stock as entry allowed!)
	datass$recruitment_time <- 1					# Time step of Recruitment
	# Spawning area: Area for recording SSB - either one area per stock or all areas
	datass$SSB_for_all_areas_pin <- 1					# SSB calculated over all areas: 1 = yes, 0 = No
	datass$spawning_all_areas <- 1					# SSB recorded for all areas combined (only useful if n_stocks=1), ignored if SSB_for_all_areas_pin = 0
	datass$spawning_areas <- "R1"					# Area where spawning occurs for each stock (one for each stock), ignored if SSB_for_all_areas_pin = 1
	datass$spawning_time <- 1					# Time step of spawning
	datass$spawning_part_mort <- 0.5					# Part of mortality that spawning fish undergo in this time step before SSBs are calculated (e.g. 0.5 = mid-step, 1=end of step)
	datass$spawning_p <- 1.0					# Spawning proportion to estimate SSB for each stock
	datass$spawning_use_total_B <- "False"				# SSB as total B in spawning area rather than mature B (false = mature B)
	datass$aging_time <- 1					# Time step when ageing occurs
	datass$growth_props <- c(0)					# Prop of annual growth that has occurred by start each time step, default = 0, must be 0 in time step of ageing
	datass$M_props <- c(1) 	 			# Prop of annual natM that occurs in each time step
	datass$baranov <- "False"				# Use Baranov equation (apply F and M simultaneously)
	datass$n_migrations <- 0					# Number of migrations
	## Migrations: First distribution of recruitment, then other migrations
	if (length(datass$regions) > 1) {
		datass$n_migrations <- 1					# Number of migrations
		datass$migration_names <- c("MoveRec_1_2")		# Defines the text labels of the migrations
		datass$migration_times <- c(1)					# Time step of each migration
		datass$migrate_from <- c("R1")				# Area from which each migration departs
		datass$migrate_to	<- c("R2")				# Area where each migration arrives
		datass$migrators <- c("all")				# Defines the fish that can migrate, either mature, immature, or all
		datass$rates_all[[1]] <- c("logistic_capped", 10, 5, 0.33)
		datass$migration_est_pin <- TRUE					# Should the model estimate migration?
		datass$migration_rates_all_low <- c( 1, 0.1, 0.01)		# logistic_capped
		datass$migration_rates_all_upp <- c(30,  10, 1.00)		# logistic_capped
		datass$migration_rates_all_prior <- "uniform"
		# datass$rates_all[[1]]	 			<- c("allvalues ",0.5,rep(0,(datass$age[3]-1))) 	# Proportion of applicable fish that migrate
		# datass$migration_rates_all_low	<- c(rep(0.01,length(cpara$age[1]:cpara$age[2])))	# allvalues
		# datass$migration_rates_all_upp	<- c(rep(1,length(cpara$age[1]:cpara$age[2])))		# allvalues
		# datass$migration_obs_pin			<- TRUE			# Are there migration observations to fit to?
		# datass$migration_AEM		<- "False"				# In est@proportions_migrating: Should ageing error be applied?	(ageing_error)
		# datass$migration_do_boot 	<- "True"				# In est@proportions_migrating: Print out parametric bootstrap
	}
	#### Model recruitment
	datass$years_to_recruit <- 1					# Years to exclude for last_free (bad Rec estimate due to too little data at end of period)
	datass$rec_y_enter <- 1					# Number of years after which a year class enters the partition
	datass$rec_standardise_YCS <- "True"				# Defines YCS to use the Haist parameterisation in Section 5.4.2
	datass$rec_YCS_years <- datass$years - datass$rec_y_enter				# Years for which YCS are provided
	#datass$rec_YCS <- as.vector(apply(mod$YCS,2,sum))/om$rec_mu		# Real Year class strengths for the stock
	datass$rec_YCS <- rep(1,length(datass$years)) 						# Year class strengths for the stock
	datass$rec_first_free <- datass$year[1] - datass$rec_y_enter				# Defines YCS range [first_free ? last_free] to estimate R0 (average)
	datass$rec_last_free <- datass$year[2] - datass$rec_y_enter - datass$years_to_recruit - 2	# Fish have to be at least years_to_recruit+2y old
	datass$rec_first_random_year <- datass$rec_last_free + 1 						# Defines the first unknown YCS
	datass$rec_year_range <- c(datass$rec_first_free, datass$rec_last_free)	# Year range from which randomised YCS are resampled
	datass$rec_year_range_N <- length(datass$rec_first_free:datass$rec_last_free)
	datass$rec_SR <- "BH"					# Stock-recruitment relationship, e.g. BH or Ricker
	datass$rec_sigma <- 0.6					# Standard deviation on the log scale of randomised YCS
	datass$rec_steepness <- 0.7					# Steepness parameter of the stock-recruitment relationship
	datass$rec_rho <- 0.0					# Lag-1 log-scale autocorrelation of randomised YCS
	datass$rec_p_male <- 0.5					# Proportion of recruits that are male
	datass$rec_randomisation_method <- "lognormal"			# Randomisation method for recruitment variability in stochastic simulations and projections
	#### Biological parameters
	# For by_sex = True: partition growth, WL, natM and maturity for females and males
	datass$size_at_age_type <- "von_Bert"
	datass$size_at_age_dist <- "normal"
	datass$estgrowth 	<- list()
	# Von Bertalanffy function (Linf, K, t0, cv)
	datass$estgrowth[[1]] <- c(2870.8, 0.02056, -4.2897, 0.1)
	datass$estWL <- list()
	# Weight-Length relationship (WLa, WLb) mm to kg
	datass$estWL[[1]] <- c(2.59e-12, 3.2064)
	datass$verify_size_weight <- c(500, 0.5, 1.5) 		# A fish of 500(mm) is between 0.5-1.5 kg
	datass$estnatM <- list()
	datass$estnatM[[1]] <- 0.155
	datass$estmaturity <- list()
	datass$estmaturity[[1]] <- c(11, 17)
	# "logistic" (x50, x to 95)
	# "ramp" (low 0, low 1)
	# "provide" (provide a vector with all values)
	datass$estpin.mat <- "ramp"
	if (datass$estpin.mat == "logistic") datass$maturity_props_all <- c("allvalues ",round(ogive_logistic(datass$ages,datass$estmaturity[[1]]),4))
	if (datass$estpin.mat == "ramp") datass$maturity_props_all <- c("allvalues ",round(ogive_ramp(datass$ages,datass$estmaturity[[1]]),4))
	if (datass$estpin.mat == "provide") datass$maturity_props_all <- c("allvalues ",round(ogive_provide(datass$ages,datass$estmaturity[[1]]),4))
	# datass$maturity_props_all		<- "allvalues_bounded 11 17  0.0  0.16 0.31 0.5000 0.69 0.84 1.0000"
	# maturity_props_female			<- "female allvalues_bounded 11 17  0.0  0.16 0.31 0.5000 0.69 0.84 1.0000"
	# maturity_props_male			<- "male allvalues_bounded   11 17  0.0  0.16 0.31 0.5000 0.69 0.84 1.0000"
	#											   Age range (min, max), then values for ages within the age range"
	## Movement
	##* I don't think movement is used in this example (because it is turned off)
	datass$move_rules <- om$move_rules
	#### Fishery parameters
	datass$list_fishery <- datass$Fish[,"Fishery"]			# Prepare lists
	datass$list_season <- datass$Fish[,"Season"]
	datass$list_region <- datass$Fish[,"Region"]
	datass$list_sel <- datass$Fish[,"Sel"]
	datass$list_q <- datass$Fish[,"q"]
	datass$selectivity_names <- levels(factor(datass$list_sel))	# Unique Sel
	## switch to control estimation of selectivity in CASAL estimation file
	datass$estimate_selectivity <- datass$selectivity_names
	datass$qq_names <- levels(factor(datass$list_q))	# Unique qs
	## Selectivity values: List with shape and starting parameters for estimation of all selectivity functions as in datass$selectivity_names
	# For examples: "double_normal_plateau 3 4 1 4 1", "double_normal 4 1 7.05", "logistic 9 4", "logistic_capped 9 4 0.7"
	# In population.csl: starting parameters
	datass$selN_all <- list()
	datass$selN_all[[1]] <- c("double_normal",  5, 2, 5)	# SelTrawl	3,4,1
	datass$selN_all[[2]] <- c("double_normal", 10, 3, 8)	# SelLL
	# In Estimation.csl: bounds								    # Lower bound, upper bound, prior (uniform, uniform-log, lognormal etc)
	datass$est_selN_all <- list()
	datass$est_selN_all[[1]] <- list(c(1, 0.1, 0.1), c(10, 10, 20), "uniform") 	# SelTrawl
	datass$est_selN_all[[2]] <- list(c(1, 0.1, 0.1), c(20, 10, 20), "uniform")  # SelLL
	#datass$est_selN_all <- list(); datass$est_selN_all[[1]] <- list(c(1, 0.1, 0.02, 1.00, 1.00), c(10, 20, 20, 12, 1.00), "uniform")  	# SelLL
	## Sex-based:
	# datass$selN_female <- list(); datass$selN_female[[1]]  <- c("double_normal_plateau", 3, 4, 1, 4, 1)
	# datass$selN_male   <- list(); datass$selN_male[[1]] 	 <- c("double_normal_plateau", 3, 4, 1, 4, 1)
	# Also change selectivity[xxx].all in @estimate (estimation.csl)
	## Catchability values: Vector with names for each catchability as in datass$qq_names
	#  If catchability should not be estimated (e.g for surveys), set lower bound = upper bound
	# In Estimation file: Bounds						# Lower bound, upper bound, prior (uniform, uniform-log, lognormal...), mu & cv (for lognormal)
	datass$qqvalues <- list()
	datass$qqvalues[[1]] <- c(1, 1, "uniform", 1, 0.5)
	datass$qqvalues[[2]] <- c(1, 1, "uniform", 1, 0.5)
	datass$qqvalues[[3]] <- c(1, 1, "uniform", 1, 0.5)
	## Future catches in projections (initial values, will be overwritten later)
	datass$future_constant_catches <- rep(200, length(datass$list_fishery))		# Length and order must match list_fishery
	# Umax
	datass$U_max <- 0.95
	###*** consider moving to doc file
	#### Sampling matrix: Define data types and years by fisheries for the assessment
	#### Data collection system
	# Determine years of sampling (not season specific - assumes that all seasons with a fishery are sampled!):
	# 1. catchage:  Age composition samples
	# 2. catchlen:  Length composition samples
	# 3. surveyage: Age survey
	# 4. catchlen:  Length survey
	# 5. index:     Index of abundance
	# 5. tagging:   Years with tag-releases
	#
	# datass$sample_years is set up for all years of the OM. No need to define season & region (they are fishery-specific)
	types <- c("catchage","catchlen","survey","surveyage","surveylen","cpue","tagging")
	datass$sample_years <- array(0, dim=c(length(types),length(datass$years),length(datass$list_fishery)),
							dimnames=list(types=types,year=datass$years,fishery=datass$list_fishery))
	# Choose the data that will be included = 1 (0 if not included)
	ycurr <- which(datass$years==Yr_last)
	age_y <- which(datass$years==2000):ycurr
	surv_y <- seq(which(datass$years==2000),ycurr,2)
	datass$sample_years["catchlen", (ycurr-5):ycurr,] <- 0
	datass$sample_years["catchage", age_y,] <- 1
	datass$sample_years["survey", surv_y,2] <- 0
	datass$sample_years["surveyage", surv_y,] <- 0
	datass$sample_years["surveylen", (ycurr-4):ycurr,] <- 0
	datass$sample_years["cpue", (ycurr-5):ycurr,] <- 0	# Relative index (cpue)
	datass$sample_years["tagging", (ycurr-5):(ycurr-1),] <- 1	# No tags in last year
	#### Abundance index, Survey data, Catch-at-age and catch-at-length
	for (ff in 1:length(datass$list_fishery)) {
	  ##** change this to be by fishery
	  ## also update the om fishery sampling here
		fish <- datass$list_fishery[ff]					# Assessment fishery
		## CPUE: Index of relative Abundance
		#datass$cpue_index[[fish]] <- NA			# Space holder
		#datass$cpue_cv[[fish]]				<- NA			# Space holder
		#datass$survey_index[[fish]]		<- NA			# Space holder
		#datass$survey_ind_var[[fish]]		<- NA			# Space holder
		#datass$survey_age_n[[fish]]		<- NA			# Space holder
		#datass$survey_len_n[[fish]]		<- NA			# Space holder
		#datass$survey_age_var[[fish]]		<- NA			# Space holder
		#datass$survey_len_var[[fish]]		<- NA			# Space holder
		#datass$catch_age[[fish]]			<- NA
		#datass$catch_len[[fish]]			<- NA
		datass$cpue_cvpar <- 0.234 		# Variance estimates used in assessment
		datass$cpue_cv_dist <- "lognormal"	# Error distribution of estimated abundance index in units (normal, lognormal) (used in assessment)
		datass$cpue_cv_process_error <- 0			# If 0, no process error is applied (= no down-weighting)
		## Survey abundances-at-age and length
		datass$surveyindexvar <- 0.3
		datass$surveyagevar <- 0.3
		datass$surveylenvar <- 0.0
		datass$survey_dist <- "lognormal"	# Error distribution of survey-at-age and survey-at-length data
		datass$survey_cv_process_error <- 0			# If 0, no process error is applied (= no down-weighting)
		## Catch-at-length and catch-at-age
		datass$catch_at_dist <- "multinomial" 	# Error distribution of catch-at-age and catch-at-length data
		datass$catch_at_r <- 1e-11			# Robustification parameter used in binomial likelihood (relevant when fitted prop is close to 1 or 0)
		## To alter the weighting of these data in the likelihood function without changing the quality of the data
		datass$pin_diff_N_in_assessment <- 1				# If 1, use these N for assessment, if 0, use N from sampling (in OM) in assessment
	}
	# For abundance, survey and tag-recapture observations (default 0.5)
	datass$proportion_mortality <- 0.5	# Proportion of step?s mortality, prior to when observations occur
	#### Tag release and recapture
	## Need to translate tagging data from OM to tagging data in Assessment (mainly account for potentially diff regions)
	## Tagging release: Releases by fishery (each fishery is region & season-specific)
	# Select tag releases from all fisheries & all years (even if not specified in sample.years):
	# for (ff in 1:length(datass$list_fishery)) {
	# 	fish	<- datass$list_fishery[ff]					# Assessment fishery
	# 	datass$tag_numbers[[fish]]   	<- NULL				# Space holder
	# 	datass$tag_props_all[[fish]]	<- NULL				# Space holder
	# }
	# datass$landings_n_len_sum <- NULL				# Scanned numbers = Landing numbers by length summed across all fisheries in a region
	datass$tag_shedding_rate <- 0.0084			# Tag shedding rate to apply to all tagging partitions
	datass$tag_release_type	<- "deterministic"	# Method for determining proportions-at-age release in an age-based model
	datass$tag_sex <- "both"			# Sex that the tagging event applies to
	#datass$tag_stock <- datass$stock_names		# Fish stock of tagging
	datass$tag_mature_only <- "False"			# Does the tagging event apply to mature or all fish?
	datass$tag_loss_props <- rep(1/datass$season[3],datass$season[3])	# Proportion of tag loss that has occurred by each time step
	datass$nogrowth_period <- 0.0				# Period of no growth after tagging applied to each tag partition
	datass$tag_mortality <- 0.1				# Proportion of tagged fish that are removed immediately after tagging (not a rate!)
	# In Estimation file
	datass$tag_proportion_scanned <- 1.00				# Proportion of the catch that is scanned (to calculate 'scanned_[year]' numbers)
	datass$tag_detection_probability <- 1.00				# Tag detection rate
	datass$tag_dispersion <- 1.0				# Dispersion term ? for the tag-recap likelihood (which is modified by multiplied by 1/?)
	datass$tag_do_bootstrap <- "True"			# Print parametric bootstraps (default: true)
	datass$tag_r <- 1e-11			# Robustification parameter used in binomial likelihood (relevant when fitted prop is close to 1 or 0)
	datass$tag_sampling_type <- "size"			# Sampling type: size: 	   Length for recaptures and tagged fish
	# size-age: Ages for tagged fish, length for scanned fish (Not implemented!!!!)
	datass$tag_y_liberty <- 6				# Number of years of liberty after tagging event for which recaptures are included in model
	#### Estimation file: Parameters
	## Estimation parameters
	datass$estimator <- "Bayes"
	datass$max_iters <- 1600
	datass$max_evals <- 10000
	datass$grad_tol <- 0.001
	datass$MCMC_start <- 0
	datass$MCMC_length <- 1500000
	datass$MCMC_keep <- 1000
	datass$MCMC_stepsize <- 0.01
	datass$MCMC_adaptive_stepsize <- "True"
	datass$MCMC_adapt_at <- c(100000, 200000, 300000, 400000)
	datass$MCMC_burn_in	<- 500
	## Estimation profiling
	datass$profile_parameter <- "initialization.B0"
	datass$profile_n <- 11
	datass$profile_l <- 40000		# 120000
	datass$profile_u <- 130000		# 160000
	# Estimation: General
	datass$ageing_error <- c("normal", 0.0)		# Ageing error type (e.g. normal) and cv of misclassification
	datass$q_method <- "nuisance"			# Estimation: q method
	# Include (1=yes, 0=no), lower bound, upper bound, prior shape, name
	datass$estim_initialization.B0 <- c(1, datass$initialB0[3], datass$initialB0[4], "uniform", "initialization.B0")
	# If B0 not estimated:
	#datass$estim_initialization.B0 <- c(0, datass$B0, datass$B0, "uniform", "initialization.B0")
	datass$estim_size_at_age.cv <- c(0, 0.05, 0.20, "uniform", "size_at_age.cv")
	datass$estim_natural_mortality.all <- c(0, 0.10, 0.25, "uniform","natural_mortality.all")
	datass$estim_natural_mortality.ogive_all <- c(0, 0.10, 0.25, "uniform", "natural_mortality.ogive_all")
	datass$estim_recruitment.YCS <- list()
	datass$estim_recruitment.YCS[[1]] <- 1		# Should YCS be estimated? 1 = Yes, 0 = No
	datass$estim_recruitment.YCS[[2]] <- c(rep(0.001, datass$rec_year_range_N),rep(1,length(datass$rec_YCS)-datass$rec_year_range_N))
	datass$estim_recruitment.YCS[[3]] <- c(rep(  100, datass$rec_year_range_N),rep(1,length(datass$rec_YCS)-datass$rec_year_range_N))
	datass$estim_recruitment.YCS[[4]] <- "lognormal"
	datass$estim_recruitment.YCS[[5]] <- "recruitment.YCS"
	datass$estim_recruitment.YCS[[6]] <- rep(1, length(datass$rec_YCS))			# mu
	datass$estim_recruitment.YCS[[7]] <- rep(1.0, length(datass$rec_YCS))			# cv
	# Estimation boundaries for Catchability and Selectivity: see above
	# Estimate: Process error					# Include (1=yes, 0=no), lower bound, upper bound, prior shape
	datass$estim_abund_cv_process_error	<- c(0, 0.001, 5.00, "uniform")
	datass$estim_survA_cv_process_error	<- c(0, 0.001, 5.00, "uniform")
	datass$estim_survS_cv_process_error	<- c(0, 0.001, 5.00, "uniform")
	# Penalty
	datass$catch_limit_penalty_multiplier <- 1000
	datass$catch_limit_penalty_log_scale <- "True"
	datass$fish_tagged_penalty_multiplier	<- 1
	datass$vector_average_penalty <- c(1, "meanYCS_1", "recruitment.YCS", 1, 100)	# Include (yes=1, No=0), label, vector, average (k), multiplier
	#### Output file: Parameters
	# Print: Estimation section
	datass$output[["print"]]$parameters <- "true"
	datass$output[["print"]]$covariance <- "true"
	datass$output[["print"]]$fits_every_eval <- "false"
	datass$output[["print"]]$objective_every_eval <- "false"
	datass$output[["print"]]$parameters_every_eval <- "false"
	datass$output[["print"]]$parameter_vector_every_eval <- "false"
	datass$output[["print"]]$fits <- "true"
	datass$output[["print"]]$resids <- "true"
	datass$output[["print"]]$pearson_resids <- "true"
	datass$output[["print"]]$normalised_resids <- "true"
	datass$output[["print"]]$estimation_section <- "true"
	# Print: Population section
	datass$output[["print"]]$requests <- "false"
	datass$output[["print"]]$initial_state <- "false"
	datass$output[["print"]]$state_annually <- "false"
	datass$output[["print"]]$state_every_step <- "false"
	datass$output[["print"]]$final_state <- "false"
	datass$output[["print"]]$results <- "false"
	# Print: Output section
	datass$output[["print"]]$yields <- "true"
	datass$output[["print"]]$unused_parameters <- "true"
	# Quantities
	datass$output[["quantities"]]$all_free_parameters <- "true"
	datass$output[["quantities"]]$nuisance_qs <- "true"
	datass$output[["quantities"]]$B0 <- "true"
	datass$output[["quantities"]]$R0 <- "true"
	datass$output[["quantities"]]$SSBs <- "true"
	datass$output[["quantities"]]$YCS <- "true"
	datass$output[["quantities"]]$true_YCS <- "true"
	datass$output[["quantities"]]$recruitments <- "true"
	datass$output[["quantities"]]$fishing_pressures <- "true"
	datass$output[["quantities"]]$actual_catches <- "false"
	#datass$output[["quantities"]]$ogive_parameters <- selectivity[Sel_Survgrp1].all

	for (ee in 1:length(datass$regions)) {
		Natage	<- paste("numbers_at[Numbers_at_age_",datass$regions[ee],"]",sep="")
		datass$output[[Natage]]$command <- "numbers_at"
		datass$output[[Natage]]$value <- paste("Numbers_at_age_",datass$regions[ee],sep="")
		datass$output[[Natage]]$step <- datass$Fish[1,"Season"]     # Season: Take Season from first fishery
		datass$output[[Natage]]$area <- datass$regions[ee]
		datass$output[[Natage]]$proportion_mortality <- datass$proportion_mortality
		datass$output[[Natage]]$years <- datass$years
	}
	# Number of projections (applies only to projections from point estimate, not from sample of posterior)
	datass$output[["n_projections"]] <- 100
	########################
				#@abundance total_biomass
				# output quantity: total biomass in all areas
				#biomass true
				#all_areas true
				#step 2
				#proportion_mortality 0.5
				#years 1970 1971 1972 1973 1974 1975 1976 1977 ? 1998 1999 2000
	      #
				#@abundance FishingVulnerableBiomass
				#biomass true
				#step 2
				#proportion_mortality 0.5
				#ogive FishingSel
				#years 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998
				#1999 2000 2001 2002 2003 2004 2005 2006 2007 2008
	#########################
	## Return para
	para[["ass"]] <- datass
	return(para)
}

#' Get CASAL data
#'
#' Get CASAL data
#' @param datass Assessment parameters and data
#' @param Yr_current Current year
#' @param om Operating model parameters
#' @param sampling Sampling parameters
#' @param obs Observations
#' @param tag tagging data
#' @param mod Results from operating model
#' @export
get_casal_data <- function(datass, Yr_current, om, sampling, obs, tag, mod) {
	## Update CASAL assessment parameters  and add data to para[["ass"]]
	# datass	<- para[["ass"]]
	# om		<- para[["om"]]
	# sampling	<- para[["sampling"]]
	# obs		<- res[["obs"]]
	# tag		<- res[["tag"]]
	# mod		<- res[["mod"]]
	datass$Yr_current <- Yr_current		# The current/assessment year (eg. 2010)
	#### Update years
	datass$year <- c(om$year[1],Yr_current,which(om$years==Yr_current))	# first, current, which year in om
	datass$years <- seq(datass$year[1],datass$year[2])					# Sequence of all years
	yrs <- 1:datass$year[3]										# Used to define years in some arrays below
	ycurr <- which(datass$years==datass$Yr_current)				# Number of current year
	## Update recruitment years
	datass$rec_YCS_years <- datass$years - datass$rec_y_enter				# Years for which YCS are provided
	#datass$rec_YCS <- as.vector(apply(mod$YCS,2,sum))/om$rec_mu		# Real Year class strengths for the stock
	datass$rec_YCS <- rep(1,length(datass$years)) 						# Year class strengths for the stock
	datass$rec_first_free <- datass$year[1] - datass$rec_y_enter				# Defines YCS range [first_free ? last_free] to estimate R0 (average)
	datass$rec_last_free <- datass$year[2] - datass$rec_y_enter - datass$years_to_recruit - 2	# Fish have to be at least years_to_recruit+2y old
	datass$rec_first_random_year <- datass$rec_last_free + 1 						# Defines the first unknown YCS
	datass$rec_year_range <- c(datass$rec_first_free, datass$rec_last_free)	# Year range from which randomised YCS are resampled
	datass$rec_year_range_N <- length(datass$rec_first_free:datass$rec_last_free)
	#### Add catches (match catches of om and assessment)
	## Observed Catch for assessment: (used in population & estimation file)
	for (ff in 1:length(datass$list_fishery)) {
		fish <- datass$list_fishery[ff]			# Assessment fishery
		seas <- datass$Fish[ff,"Season"]
		omfish <- unique(datass$match_fishery["OM",] [datass$match_fishery["Ass",] %in% datass$Fish[ff,"Fishery"]])
		omreg	<- unique(datass$match_region["OM",] [datass$match_region["Ass",] %in% datass$Fish[ff,"Region"]])
		# Observation error is added in OM (hence obs_catch is used here rather than mod$real_catch)
		# Data from first OM fishery matching this assessment fishery
		datass$catch[[fish]] <- apply(mod$obs_catch[[omfish[1]]][,yrs,,seas,omreg,drop=FALSE],c(1,2,3,4),sum)
		if(length(omfish) >1)		# If there are more than one OM fisheries matching this assessment fishery
			for (i in 2:length(omfish))
				datass$catch[[fish]] <- datass$catch[[fish]] + apply(mod$obs_catch[[omfish[i]]][,yrs,,seas,omreg,drop=FALSE],c(1,2,3,4),sum)
	}
	## Checks:
	# All fisheries are represented
	if(length(om$fishery) > length(datass$list_metier))
		print(paste("Metiers in OM that are not included in Assessment: ",om$fishery[! om$fishery %in% datass$list_metier],sep=""))
	# print("")
	# print("Fisheries in OM and Assessment:")
	# print(datass$Fish)
	# print("")
	# print("Real catch observed in OM by fishery, season and region: ")
	# print(apply(om$catch,2:4,sum))
	# print("Real catch observed in OM by fishery, season and region: ")
	# print(apply(om$catch,2:4,sum))
	# print("Real catch observed in Assessment by fishery, season and region: ")
	# print(apply(om$catch,2:4,sum))
	#### Abundance index, Survey data, Catch-at-age and catch-at-length
	for (ff in 1:length(datass$list_fishery)) {
		fish <- datass$list_fishery[ff]		# Assessment fishery
		seas <- datass$Fish[ff,"Season"]
		omfish <- unique(datass$match_fishery["OM",] [datass$match_fishery["Ass",] %in% datass$Fish[ff,"Fishery"]])
		omreg <- unique(datass$match_region["OM",] [datass$match_region["Ass",] %in% datass$Fish[ff,"Region"]])
		# Copy data from first match OM Fishery and region - omfish[1] and omreg[1]
		datass$cpue_index[[fish]] <- obs[[omfish[1]]]$cpue[,yrs,,seas,omreg[1],drop=FALSE]
		datass$survey_index[[fish]] <- obs[[omfish[1]]]$survey_index[,yrs,,seas,omreg[1],drop=FALSE]
		datass$survey_age_n[[fish]] <- obs[[omfish[1]]]$survey_age_n[,yrs,,seas,omreg[1],drop=FALSE]
		datass$survey_len_n[[fish]] <- obs[[omfish[1]]]$survey_len_n[,yrs,,seas,omreg[1],drop=FALSE]
		datass$catch_age[[fish]] <- obs[[omfish[1]]]$age_sample_n[,yrs,,seas,omreg[1],drop=FALSE]
		datass$catch_len[[fish]] <- obs[[omfish[1]]]$len_sample_n[,yrs,,seas,omreg[1],drop=FALSE]
		# For catch weighting: Use total catch per year (pooled over sexes)
		fish_catch <- apply(datass$catch[[fish]],2,sum)
		# If there are more than one OM fishery matching this assessment fishery, add observations from the other OM fisheries
		if(length(omfish) > 1) {
			# Adjust/Catch-weight age and length data for omfish[1] (survey and cpue are effectively catch-weighted already)
			catchwt <- sweep(mod$obs_catch[[omfish[1]]][,yrs,,seas,omreg[1],drop=FALSE],c(1:3),fish_catch,"/")
			datass$cpue_index[[fish]] <- sweep(datass$cpue_index[[fish]],2:5,catchwt,"*")	# Adjust/Catch-weight data from omfish[1] first
			datass$survey_index[[fish]] <- sweep(datass$survey_index[[fish]],2:5,catchwt,"*")
			# datass$survey_age_n[[fish]]	<- datass$survey_age_n[[fish]]		# No adjustment needed
			# datass$survey_len_n[[fish]]	<- datass$survey_len_n[[fish]]
			datass$catch_age[[fish]] <- sweep(datass$catch_age[[fish]],2:5,catchwt,"*")
			datass$catch_len[[fish]] <- sweep(datass$catch_len[[fish]],2:5,catchwt,"*")
			for (i in 2:length(omfish)) {
				catchwt <- sweep(mod$obs_catch[[omfish[i]]][,yrs,,seas,omreg[i],drop=FALSE],c(1:3),fish_catch,"/")
				datass$cpue_index[[fish]] <- datass$cpue_index[[fish]]   + sweep(obs[[omfish[i]]]$cpue[,yrs,,seas,omreg[i],drop=FALSE],2:5,catchwt,"*")
				datass$survey_index[[fish]] <- datass$survey_index[[fish]] + sweep(obs[[omfish[i]]]$survey_index[,yrs,,seas,omreg[i],drop=FALSE],2:5,catchwt,"*")
				datass$survey_age_n[[fish]] <- datass$survey_age_n[[fish]] + obs[[omfish[i]]]$survey_age_n[,yrs,,seas,omreg[i],drop=FALSE]
				datass$survey_len_n[[fish]] <- datass$survey_len_n[[fish]] + obs[[omfish[i]]]$survey_len_n[,yrs,,seas,omreg[i],drop=FALSE]
				datass$catch_age[[fish]] <- datass$catch_age[[fish]]    + sweep(obs[[omfish[i]]]$age_sample_n[,yrs,,seas,omreg[i],drop=FALSE],2:5,catchwt,"*")
				datass$catch_len[[fish]] <- datass$catch_len[[fish]]    + sweep(obs[[omfish[i]]]$len_sample_n[,yrs,,seas,omreg[i],drop=FALSE],2:5,catchwt,"*")
			}
		}
		## CPUE: Index of relative Abundance
		datass$cpue_cv[[fish]] <- datass$cpue_index[[fish]] 	# Create array
		datass$cpue_cv[[fish]][] <- datass$cpue_cvpar 			# Variance estimates used in assessment
		## Survey abundances-at-age and length
		datass$survey_index_var[[fish]] <- datass$survey_index[[fish]]	# Create array
		datass$survey_age_var[[fish]] <- datass$survey_age_n[[fish]] 	# Create array
		datass$survey_len_var[[fish]] <- datass$survey_len_n[[fish]] 	# Create array
		datass$survey_index_var[[fish]][] <- datass$surveyindexvar
		datass$survey_age_var[[fish]][] <- datass$surveyagevar
		datass$survey_len_var[[fish]][] <- datass$surveylenvar
	}
	#### Tag release
	## Tag-release events (in assessment)
	ttag <- datass$sample_years["tagging",1:ycurr,,drop=FALSE]
	tag_names <- NULL
	for (yy in 1:length(dimnames(ttag)$year)) {
		for (ff in 1:length(dimnames(ttag)$fishery)) {
			tyear <- dimnames(ttag)$year[yy]
			tfish <- dimnames(ttag)$fishery[ff]
			treg <- as.vector(datass$Fish[datass$Fish[,"Fishery"]==tfish,"Region"])
			if(ttag[,yy,ff]>0) {
			  tag_names <- rbind(tag_names,
			                     cbind(TagName=paste("Tags",tyear,"_",treg,sep=""),
			                           TagEvent=paste("Tags",tyear,"_",treg,"_",tfish,sep=""),
			                           Year=tyear,
			                           Fishery=tfish,
			                           Region=treg))
			}
		}
	}
	datass$tag_names <- tag_names
	## Tag release data: Releases by fishery (each fishery is region & season-specific)
	# Need to translate tagging data from OM to tagging data in Assessment (mainly account for potentially diff regions)
	# Select tag releases from all fisheries & all years (even if not specified in sample.years):
	for (ff in 1:length(datass$list_fishery)) {
		fish <- datass$list_fishery[ff]			# Assessment fishery
		seas <- datass$Fish[ff,"Season"]
		omfish <- unique(datass$match_fishery["OM",] [datass$match_fishery["Ass",] %in% datass$Fish[ff,"Fishery"]])
		omreg <- unique(datass$match_region["OM",] [datass$match_region["Ass",] %in% datass$Fish[ff,"Region"]])
		# The sum of length and age tag_N may not be equal due to rounding when converting ages to lengths.
		# Since the assessment model uses length data, length data is used to extract total numbers tagged
		# Copy data from first match OM Fishery and region - omfish[1] and omreg[1]
		tagNumLen <- obs[[omfish[1]]]$tag_len_n[,,,seas,omreg[1],drop=FALSE]		# Sum by assess region
		tagNum <- apply(obs[[omfish[1]]]$tag_len_n[,,,seas,omreg[1],drop=FALSE],2,sum)	# By year
		# If there are more than one OM fishery matching this assessment fishery, add observations from the other OM fisheries
		if(length(omfish) > 1) {
			for (i in 2:length(omfish)) {
				tagNumLen <- tagNumLen + obs[[omfish[i]]]$tag_len_n[,,,seas,omreg[i],drop=FALSE]		# Sum by assess region
				tagNum <- tagNum + apply(obs[[omfish[i]]]$tag_len_n[,,,seas,omreg[i],drop=FALSE],2,sum)			# By year
			}
		}
		# Store data and convert to proportions
		datass$tag_numbers[[fish]] <- tagNum
		datass$tag_props_all[[fish]] <- sweep(tagNumLen,c(2),tagNum,"/")
		datass$tag_props_all[[fish]][is.na(datass$tag_props_all[[fish]])] <- 0
		# Check: apply(datass$tag_props_all[[fish]],2,sum)
	}
	## Tag Recapture data: Recapture by Release area * recapture area
	# Prepare storage
	datass$Rec_events <- NULL
	datass$Rec_data	<- list()
	datass$scannedN	<- list()
	# Prepare growth for tag$iTags
	ss <- as.numeric(as.character(datass$Fish[1,"Season"])) 	# Season: Take Season from first fishery
	ssages <- om$ages + om$growth_props[ss]
	# Convert tag_names to unique release areas
	tag_names	<- tag_names[!duplicated(tag_names[,"TagName"]),c("TagName","Year","Region")]
	for (ff in 1:nrow(tag_names)) {
		# Retrieve Tagging year and area
		tyear	<- as.numeric(tag_names[ff,"Year"])
		tyy	<- which(datass$years  %in% tyear)			# Year number
		treg <- tag_names[ff,"Region"]
		tr_om <- unique(datass$match_region["OM",] [datass$match_region["Ass",] %in% treg])
		# Define years of recaptures: (tyear+1):(tyear+datass$tag_y_liberty)
		ryear <- (tyear+1) : min(tyear+datass$tag_y_liberty,(datass$year[2]))
		ryy <- which(datass$years %in% ryear)			# Year number
		# For each recapture area in the assessment:
		for (rr in 1:length(datass$regions)){					# rr = Assessment Recapture regions
			# Determine equivalent om regions
			rr_om	<- unique(datass$match_region["OM",] [datass$match_region["Ass",] %in% datass$regions[rr]])	# Recapture region OM
			## For Pooled tags
			if (sum(tag$recaps_len) > 0) {	# If any data in recaps_len
				Tpool	<- tag$recaps_len[,ryy,,,rr_om,tyy,tr_om,drop=FALSE]
				# Sum recapture Numbers over by length, year, and sex
				dat <- round(apply(Tpool,drop=FALSE,c(1,2,3),sum),1)
				if(datass$by_sex == 0) dat	<- apply(dat,c(1,2),sum)	# Sum over sex
			}
			## For Individual tags
			if (nrow(tag$iTags) > 1) {	# If any data in tag$iTags
				Tind <- tag$iTags[tag$iTags$RelY %in% tyear & tag$iTags$RelArea %in% tr_om & 	# Select tags with F = 1
									 tag$iTags$LastY %in% ryear & tag$iTags$LastArea %in% rr_om & tag$iTags$F == 1,]
				# Sum recapture Numbers over by age, year, and sex
				dat_age <- xtabs(~LastAge + LastY + Sex, data=Tind)
				dat <- array(0,dim=c(length(sampling$len_classes),length(ryear),length(om$sex)),	# Create length-specific array
								   dimnames=list(len=sampling$len_classes,year=ryear,sex=om$sex))
				if(sum(dat_age) >0) {	# If 0, then dat will contain only 0
					for (i in 1:length(om$sex)) {
						if(om$sex[i] %in% dimnames(dat_age)$Sex){
							# Fill dat_age into matrix with all ages (dat_age1)
							dat_age1	<- matrix(0,nrow=length(om$ages),ncol=length(ryear),dimnames=list(om$ages,ryear))
							for (x in 1:nrow(dat_age))
								for (y in 1:ncol(dat_age)) dat_age1[rownames(dat_age1) %in% rownames(dat_age)[x],
																	colnames(dat_age1) %in% colnames(dat_age)[y]] <- dat_age[x,y,om$sex[i]]
							## Convert to lengths for each sex (growth may be sex-specific)
							for (y in 1:ncol(dat_age1)) {
								dat[,y,i] <- round(sample_lengths(ssages,sampling$len_classes,dat_age1[,y],om$growth[[om$sex[i]]]),1)
							}
						}
					}
				}
				if(datass$by_sex == 0) dat	<- apply(dat,c(1,2),sum)
			}
			# Scanned numbers = Landing numbers by length summed across all fisheries in a region
			scannedN		 <- round(apply(mod$landings_n_len_sum[,ryy,,,rr_om,drop=FALSE],c(1,2,3),sum) * datass$tag_proportion_scanned,0)
			if (datass$by_sex == 0) scannedN <- apply(scannedN,c(1,2),sum) 	# Sum over sex
			dat[scannedN==0] <- 0	# Turn prop of tagged fish to 0 for age or length classes in which 0 fish have been scanned
			# Store name of recapture event and data
		  revent <- paste(tag_names[ff,"TagName"],"_",datass$regions[rr],sep="")
			datass$Rec_events <- rbind(datass$Rec_events,c(tag_names[ff,"TagName"],RecName=revent,RecRegion=datass$regions[rr]))
			datass$Rec_data[[revent]] <- dat
			datass$scannedN[[revent]] <- scannedN
		}
	}
	#### Estimation file: Estimated Parameters
	datass$estim_recruitment.YCS[[2]] <- c(rep(0.001, datass$rec_year_range_N),rep(1,length(datass$rec_YCS)-datass$rec_year_range_N))
	datass$estim_recruitment.YCS[[3]] <- c(rep(  100, datass$rec_year_range_N),rep(1,length(datass$rec_YCS)-datass$rec_year_range_N))
	datass$estim_recruitment.YCS[[6]] <- rep(1, length(datass$rec_YCS))			# mu
	datass$estim_recruitment.YCS[[7]] <- rep(1.0, length(datass$rec_YCS))			# cv
	#### Output file: Parameters
	for (ee in 1:length(datass$regions)) {
		Natage	<- paste("numbers_at[Numbers_at_age_",datass$regions[ee],"]",sep="")
		datass$output[[Natage]]$years <- datass$years
	}
	return(datass)
}

#' Update CASAL file names
#'
#' Update CASAL file names
#' @param params Should we call this be para for consistency?
#' @export
update_casal_file_names	<- function(params) {
  ## Updates casal file names with input/outputprefix
  inputprefix <- params$inputprefix
  outputprefix <- params$outputprefix
  # specify casal files names
  params$pop_skel_csl <- paste(inputprefix, "population_skel.csl",sep="")
  params$pop_csl <- paste(inputprefix, "population.csl",sep="")
  params$est_skel_csl <- paste(inputprefix, "estimation_skel.csl",sep="")
  params$est_csl <- paste(inputprefix, "estimation.csl",sep="")
  params$output_skel_csl <- paste(inputprefix, "output_skel.csl",sep="")
  params$output_csl <- paste(inputprefix, "output.csl",sep="")
  params$mpd_dat <- paste(outputprefix,"mpd.dat",sep="")
  params$output_log <- paste(outputprefix,"output.log",sep="")
  params$profile_dat <- paste(outputprefix,"profile.dat",sep="")
  params$profile_out <- paste(outputprefix,"profile.out",sep="")
  params$profile_log <- paste(outputprefix,"profile.log",sep="")
  params$proj_dat <- paste(outputprefix,"proj.dat",sep="")
  params$mvnsamples_dat <- paste(outputprefix,"mvn_samples.dat",sep="")
  ## return the parameters
  return(params)
}
