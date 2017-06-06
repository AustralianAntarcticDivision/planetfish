## Operating model specification
## This file contains functions that run the operating model

#' Setup operating model objects
#'
#' Setup operating model objects ('res' with pop, fleet and obs).
#' @inheritParams om
#' @export
setup_om_objects <- function(para) {
  ## create an empty list to store the results
  res <- list()
	## Extract parameters
	om <- para[["om"]]
	sampling <- para[["sampling"]]
	## Create an object to hold the biological information
	res[["pop"]] <- list(n = create_array(om, sampling, type="age"), # iniAgeQuant: Numbers at age
	                     move = create_array(om, sampling, type="move"), ## movement array
	                     fec = create_array(om, sampling, type="age"), # iniAgeQuant: Fecundity
	                     wt = create_array(om, sampling, type="age"), # iniAgeQuant: Weight
	                     m = create_array(om, sampling, type="age")) # iniAgeQuant: Natural mortality
	## Create an object to hold the fishery information
	res[["fleet"]] 	<- list()
	for (ii in 1:length(om$fishery)) #simple
		res[["fleet"]][[om$fishery[ii]]] <- list(effort = create_array(om, sampling, type="simple"), # was iniQuant
		                                         landings = create_array(om, sampling, type="simple"), # was iniQuant
		                                         landings_n = create_array(om, sampling, type="age"), # was iniAgeQuant
		                                         landings_sel = create_array(om, sampling, type="age")) # was iniAgeQuant
	## Create an object to hold the population observations
	res[["obs"]] <- list()
	for (ii in 1:length(om$fishery))
		res[["obs"]][[om$fishery[ii]]] <- list(age_sample_N = create_array(om, sampling, type="simple"), # was iniQuant # Age observations: Total number of observations
		                                       age_sample_n = create_array(om, sampling, type="age"), # was iniAgeQuant: Age observations: Numbers-at-age
		                                       len_sample_N = create_array(om, sampling, type="simple"), # was iniQuant: Length observations: Total number of observations
		                                       len_sample_n = create_array(om, sampling, type="length"), # was iniLenQuant: Length observations: Numbers-at-length
		                                       survey_index = create_array(om, sampling, type="simple"), # was iniQuant: Survey: Absolute Index of abundance (scaled to vulnerable B by sampling$survey_q)
		                                       survey_age_n = create_array(om, sampling, type="age"), # was iniAgeQuant: Survey: Numbers-at-age (scaled to vulnerable population ages by sampling$survey_q)
		                                       survey_len_n = create_array(om, sampling, type="length"), # was iniLenQuant: Survey: Numbers-at-length  (scaled to vulnerable population lengths by sampling$survey_q)
		                                       cpue = create_array(om, sampling, type="simple"), # was iniQuant		# Catch rates: Relative index of abundance (scaled to vulnerable B by sampling$cpue_q)
		                                       tag_N = create_array(om, sampling, type="simple"), # was iniQuant		# Tagged fish: Total number
		                                       tag_sel = create_array(om, sampling, type="age"), # was iniAgeQuant: Tagged fish: Tagging selectivity
		                                       tag_age_n = create_array(om, sampling, type="age"), # was iniAgeQuant: Tagged fish: Numbers-at-age of fish
		                                       tag_len_n = create_array(om, sampling, type="length")) # was iniLenQuant: Tagged fish: Numbers-at-length of fish
	## Create an object to hold the model information
	res[["mod"]] <- list(rec = create_array(om, sampling, type="simple"), # was iniQuant: Annual Recruitment at age 1
	                     ##* create these arrays using 'create_array'
	                     ##* this is a standard array but with year names - 1
	                     YCS = array(0, dim=c(1,om$n_years,length(om$sex),om$n_seasons,om$n_regions),	# Year Class Strength at age = 0
												      dimnames=list(age=0,year=om$years-1,sex=om$sex,season=om$seasons,area=om$regions)),
											 ##* this is a new array but with year names
											 init_age = array(0, dim=c(om$n_ages,1,length(om$sex),1,om$n_regions),	# Initial age structure at start of year=1 and season=1
												                dimnames = list(ages=om$ages,year=om$years[1],sex=om$sex,season=1,area=om$regions)),
											 prior_rec = 0, # Annual recruitment of year classes that contribute to initial age structure
											 totB = create_array(om, sampling, type="simple"), # was iniQuant: Mid-year total biomass
											 ssb = create_array(om, sampling, type="simple"), # was iniQuant: Mid-year SSB
											 ssb0 = 0,					# SSB0 estimate used in model
											 ssb0_casal	= 0,					# SSB0 estimate: As in CASAL model
											 ssb0_init = 0,					# SSB0 estimate: Calculated from initial age composition
											 ssb0_stoch = 0,					# SSB0 estimate: Median of stochastic realisations of SSB0 (used to be ssb0_median)
											 # Sum of fishery-specific H by age
											 h_species = create_array(om, sampling, type="age"),
											 h_fishery = list(),
											 h_fishery_sum = create_array(om, sampling, type="simple"), # was iniQuant: Sum of harvest rates (summed over fishery): sum(h_fishery) <= h_max
											 h_max = 0,					# Maximum H allowed in a season and area
											 landings_n_len	= list(),				# Landing numbers by length and fishery
											 landings_n_len_sum = create_array(om, sampling, type="length"), # was iniLenQuant: Landing numbers by length (summed over fishery)
											 real_catch = list(),				# Actual catch, by sex; similar to om$catch, but by sex
											 obs_catch = list(),				# Will be observed (or reported) catch taken
											 ##* use 'create_array'
											 ##* this is a new array, should it have a fishery dimension
											 TAC = array(NA, dim=c(1,om$n_years+1),dimnames=list(tac=1,year=c(om$years,om$n_years+1))) )	# TAC recommendation - Model output
  # loop over fisheries
	for (ff in 1:length(om$fishery)) {
		res[["mod"]]$h_fishery[[om$fishery[ff]]] <- create_array(om, sampling, type="simple") # was iniQuant
		res[["mod"]]$real_catch[[om$fishery[ff]]] <- create_array(om, sampling, type="simple") # was iniQuant
		res[["mod"]]$obs_catch[[om$fishery[ff]]] <- create_array(om, sampling, type="simple") # was iniQuant
		res[["mod"]]$landings_n_len[[om$fishery[ff]]] <- create_array(om, sampling, type="length") # was iniLenQuant:
	}
	## Create an object to hold the tagging information
	res[["tag"]] <- list(tags = create_array(om, sampling, type="tag_age"), # was iniAgeQuantT: N-at-age of fish tagged in an area available to the entire fishery
	                     recaps = create_array(om, sampling, type="tag_age"), # was iniAgeQuantT: N-at-age of fish recaptured in an area by the entire fishery
	                     recaps_len	= create_array(om, sampling, type="tag_length"), # was iniLenQuantT: N-at-len of fish recaptured in an area by the entire fishery
	                     ##* use 'create_array' or perhaps not for this array
	                     # If sampling$pin_tag_Method = "Pool": tags, recaps and recaps_len is used
	                     # tags & recaps are specified by age/len, sex, season, recap year & area, tag year & area,
	                     # recaps/recaps_len are summed over all fisheries in a year
	                     # If sampling$pin_tag_Method = "Ind": iTags is used
	                     # iTags follows individual tags: Number, sex, Release year, age, season and area (could be changed to lat and lon);
	                     # 	last year, age, season area; M = dead by natM (0/1), F = dead by fishing/recapture (0/1);
	                     #   if fish is dead by fishing, last year/age/season/area becomes recapture year/age/season/area
	                     iTags = array(0, dim=c(1,13),dimnames=list(c(""),c("Sex","RelY","RelAge","RelSeas","RelArea",
																			  "LastY","LastAge","LastSeas","LastArea","Dead","RelM","M","F")) ))	# Info on individual Tags (release year & area, current year & area, status (dead, fished)...
	return(res)
}

#' Populate om object with data
#'
#' Populate om object with data
#' @inheritParams om
#' @export
populate_om_objects <- function(para, res) {
	## Extract parameters
	om <- para[["om"]]
	sampling <- para[["sampling"]]
	pop <- res[["pop"]]
	fleet <- res[["fleet"]]
	obs <- res[["obs"]]
	mod <- res[["mod"]]
	#### Biological information (pop)
	for (ss in 1:om$n_seasons) {
		# Seasonal natural mortality
		pop$m[,,,ss,] <- om$natM * om$natM_props[ss]
		# Length and weight calculations
		ssages <- om$ages + om$growth_props[ss]
		for (i in 1:length(om$sex)) {
			pop$wt[,,i,ss,] <- calc_VBweight(ssages,pop$wt[,,i,ss,],om$growth[[om$sex[i]]],om$WL[[om$sex[i]]])   # tonnes
			# Ogives are by season and ages (process is age-based thus do not use mid-year ages)
			# if (om$pin_mat == "logistic") pop$fec[] <- ogive_logistic(om$ages,om$maturity[[om$sex[i]]])
			# if (om$pin_mat == "ramp") pop$fec[] <- ogive_ramp(om$ages,om$maturity[[om$sex[i]]])
			# if (om$pin_mat == "provide") pop$fec[] <- ogive_provide(om$ages,om$maturity[[om$sex[i]]])
			pop$fec[] <- ogive(om$pin_mat, om$ages, om$maturity[[om$sex[i]]])
		}
	}
	#### Fishery data: effort and selectivity
	for (ff in 1:length(om$fishery)) {
		# Effort split by sex (For calculations of cpue etc, total effort should be used)
		for (rr in 1:om$n_regions)
			fleet[[ff]]$effort[1,,,1:om$n_seasons,rr] <- om$effort[,ff,1:om$n_seasons,rr]
		fleet[[ff]]$effort[1,,,1:om$n_seasons,] <- sweep(fleet[[ff]]$effort[1,,,1:om$n_seasons,,drop=FALSE],3,om$eff_sex,"*")
		for (ss in 1:om$n_seasons) {rr
			# Ogives are by season, mid-year ages for certain ogives
			# Age-specific processes occur at ages 1,2,3.. (not 1.5, 2.5...)
		  fleet[[ff]]$landings_sel[] <- ogive(om$pin_sel[[ff]], om$ages, om$select[[ff]])
		}
	}
	### Observations: May have to be changed if sampling is not uniform e.g. across seasons or regions
	for (ff in 1:length(om$fishery)) {
		## Catch-at observations
		obs[[ff]]$age_sample_N[] <- sampling$catchage_N/length(om$sex)	# Total sample size for catch age
		obs[[ff]]$len_sample_N[] <- sampling$catchlen_N/length(om$sex)	# Total sample size for catch length
		## Survey abundance - no input sample sizes needed
		## Tagging
		# 1: Fixed number of tags released
		if (sampling$pin_tag_N == "Fixed") obs[[ff]]$tag_N[] <- sampling$tag_N[ff]		# Number of tags released by fishery
		# 2: Tagging rate depending on catch
		if (sampling$pin_tag_N == "Catch") obs[[ff]]$tag_N[] <- 0	# Number of tags released depends on catch, thus defined later
		## specify the tagging selectivity
		obs[[ff]]$tag_sel[] <- ogive(sampling$pin_tag_sel[[om$fishery[[ff]]]],
		                             om$ages, sampling$tag_select[[om$fishery[[ff]]]])
	}
	### Model estimates
	# Maximum harvest rate
	mod$h_max	<- om$max_harvest
	## overwrite the pop, fleet, obs and mod objects, retain other components
	res[["pop"]] <- pop
	res[["fleet"]] <- fleet
	res[["obs"]] <- obs
	res[["mod"]] <- mod
	## return res
	res
}

#' Get initial population numbers and calculate SSB0
#'
#' Get initial population numbers and calculate SSB0
#' @inheritParams om
#' @export
get_initial_pop <- function(para, res) {
	## Mid-year SSB: Apply mid-year (in season 2) weight, fecundity and 0.5*natM
	## Assuming equal numbers of all sex
	## Extract parameters
	om <- para[["om"]]
	pop <- res[["pop"]]
	mod <- res[["mod"]]
	## Initial age composition
	##* need a switch here
	switch(om$pin_init_age,
	       stoch = {
	         init_age	<- initial_ages(method="stoch", mu=om$rec_mu, sigma=om$rec_sigma,
	                                  age=om$age,M=om$natM)
	         },
	       provide = {
	         init_age <- initial_ages(method="provide", age=om$age, M=om$natM,
	                                          age_series=om$init_age_series)
	         },
	       user = init_age <- initial_ages(method="user", age=om$age,
	                                       age_series=om$init_age_series))
	#init_age <- init_age_comp(om$rec_mu, om$rec_sigma, om$age, om$natM)
	# !! Use pop$m instead of om$natM if M varies between years or regions
	## Distribute initial population numbers to sex and regions
	for (rr in 1:om$n_regions)
		for (i in 1:length(om$sex))
		  mod$init_age[,1,i,1,rr] <- init_age * om$rec_area[rr] * om$rec_sex[i]
	## Initiation of age composition to all regions
	if(om$n_regions > 1) 	mod$init_age <- move_init(n=mod$init_age, om=om)
	## Retrieve init_age prior to start of the fishery
	init_age <- apply(mod$init_age,1,sum) # Retrieve init_age (may be different from initial init_age if move_init was used)
	init_rec_prior <- init_age[(om$age[1]+1):(om$n_ages-1)] * exp(cumsum(om$natM[om$age[1]:(om$n_ages-2)]))
	mod$prior_rec	<- as.vector(rev(init_rec_prior)) # Recruitment prior to the fishery
	#### Calculate mid-year SSB0:
	midyfecwt <- pop$wt * pop$fec * exp(-0.5*pop$m)	# Fecundity * weight * mid-year M ##* rename this
	## Option "init_age": Calculate SSB0 from initial age structure at start of first year (numbers by sex & area)
	# Use actual init_age from above
	##* suggest replacing this with a function
	est_ssb_init 	<- 0
	for (rr in 1:om$n_regions)
		for (i in 1:length(om$sex))
		  est_ssb_init <- est_ssb_init + sum(init_age * om$rec_area[rr] * om$rec_sex[i] * midyfecwt[,1,i,1,rr])
	mod$ssb0_init <- est_ssb_init
	## Option "stoch": Calculate SSB0 as median of stochastic realisations of SSB0
	est_ssb_stoch <- vector(mode="numeric", length=om$niter_ssb)
	for (ii in 1:om$niter_ssb) {
		init_age <- init_age_comp(om$rec_mu, om$rec_sigma, om$age, om$natM)
		for (rr in 1:om$n_regions)
			for (i in 1:length(om$sex))
			  est_ssb_stoch[ii] <- est_ssb_stoch[ii] + sum(init_age * om$rec_area[rr] * om$rec_sex[i] * midyfecwt[,1,i,1,rr])
	}
	mod$ssb0_stoch <- median(est_ssb_stoch)
	## Option 3: Calculate SSB0 as in CASAL (using R0 and natM only) = initial SSB in CASAL
	init_age <- init_age_comp(om$rec_mu, 0.0, om$age, om$natM)		# Equilibrium = no variability in R0
	est_ssb_casal <- 0
	for (rr in 1:om$n_regions)
		for (i in 1:length(om$sex))
		  est_ssb_casal <- est_ssb_casal + sum(init_age * om$rec_area[rr] * om$rec_sex[i]	* midyfecwt[,1,i,1,rr])
	mod$ssb0_casal <- est_ssb_casal
	## Define mod$ssb0
	switch(om$B0_calculation_method,
	       casal = mod$ssb0	<- mod$ssb0_casal,
	       stoch = mod$ssb0	<- mod$ssb0_stoch,
	       init_age = mod$ssb0 <- mod$ssb0_init)
	#### Return results to mod
	res[["mod"]] <- mod
	return(res)
}

#' Run Annual OM loops with/without assessment
#'
#' Run Annual Operating model loops with/without assessment
#'
#' Time sequence
#'
#' 1. Ageing
#'
#' 2. Recruitment
#'
#' 3. Maturation
#'
#' 4. Migration
#'
#' 5. Natural and fishing mortality
#'
#' 6. Sampling (Age samples, length samples)
#'
#' 7. Tag release and recaptures
#'
#' @inheritParams om
#' @param intern CASAL is run internally (default = TRUE) with no output sent
#' to the R console or it is run with output sent to R console (FALSE, use for
#' debugging)
#' @export
run_annual_om <- function(para, res, intern = TRUE) {
	## Extract parameters
	ctrl <- para[["control"]]		# para: In only
	om <- para[["om"]]				# para: In only
	sampling <- para[["sampling"]]		# para: In only
	datass <- para[["ass"]]			# para: In only
	pop <- res[["pop"]]				# res: In/out
	fleet <- res[["fleet"]]			# res: In/out
	obs	<- res[["obs"]]				# res: In/out
	mod	<- res[["mod"]]				# res: In/out
	tag <- res[["tag"]]				# res: In/out
	##* print statements slow things down
	# print start
	print("Annual loops of the operating model")
	ptm <- proc.time()		# For time measurement
	#### Calculate population numbers, apply fishery & sample catch (by region)
	for (y in 1:om$n_years) {
		print(paste("######### Year ",om$years[y]," #########",sep=""))
		## 1. Ageing & 2. Recruitment
		if(y == 1)  		# First year: Copy numbers from mod$init_age to first year
			pop$n[,y,,1,] <- mod$init_age
		if(y > 1) {			# Subsequent years
			# Recruitment
			# Recruitment with SR-relationship (BH or Ricker) is based on last-season SSB
			# Recruitment is then distributed between regions and sex (make sure SSB0 and SSB refer to the correct sexes (f or m&f))
			init <- stock_recruit(method=om$pin_rec, mu=om$rec_mu, sigma=om$rec_sigma,
			                      SSB0=mod$ssb0, SSB=sum(mod$ssb[1,y-1,,om$n_seasons,]),
			                      rec_series=om$rec_series, h=om$rec_h, year=y)
			for (rr in 1:om$n_regions) {
				for (i in 1:length(om$sex)) pop$n[1,y,i,1,rr] <- init * om$rec_area[rr] * om$rec_sex[i]
			}
			# Move population numbers from end of last year to start of this year
			pop$n[2:om$n_ages,y,,1,] <- pop$n[1:(om$n_ages-1),y-1,,om$n_seasons,]
			## this is the plus group
			pop$n[om$n_ages,y,,1,] <- pop$n[om$n_ages,y,,1,] + pop$n[om$n_ages,y-1,,om$n_seasons,]
			# Increase all ages of tagged fish available to fishing (no fish tagged in y=0). Fish not separated by season
			if (sampling$pin_tagging) {
				if(sampling$pin_tag_Method == "Pool") {
					plusgr <- tag$tags[om$n_ages,y-1,,om$n_seasons,,,]
					tag$tags[2:om$n_ages,y,,1,,,] <- tag$tags[1:om$n_ages-1,y-1,,om$n_seasons,,,]
					tag$tags[om$n_ages,y,,1,,,] <- tag$tags[om$n_ages,y,,om$n_seasons,,,] + plusgr
					tag$tags[1,y,,1,,,] <- 0
				}
				if(sampling$pin_tag_Method == "Ind") {
					liveT	<- which(tag$iTags[,"Dead"] == 0)
					tag$iTags[liveT,"LastY"] <- tag$iTags[liveT,"LastY"] + 1			# Update latest year
					tag$iTags[liveT,"LastAge"] <- tag$iTags[liveT,"LastAge"] + 1		# Update latest age
					tag$iTags[liveT,"LastSeas"]	<- 1									# Update latest season (first season)
					tag$iTags[tag$iTags[,"LastAge"] > om$n_ages,"LastAge"] <-om$n_ages	# plusgroup
				}
			}
		}
		# Store recruitment
		mod$rec[1,y,,1,] <- pop$n[1,y,,1,]
		# Back-calculate YCS from recruitment at age=1 to age=0 assuming annual natM of age=1
		for (rr in 1:om$n_regions)
			if(om$n_seasons==1) {
				totM <- pop$m[1,y,,,rr]
			} else {
				totM <- apply(pop$m[1,y,,,rr,drop=FALSE],c(4),sum)
			}
		mod$YCS[1,y,,1,rr] <- pop$n[1,y,,1,rr] * exp(totM*om$age[1])	# Store YCS
		## Catches for each year, season and region:
		## 1. Catch in om$catch is used as annual catch for each fleet
		## However, TAC is used instead if there is a TAC for the year in mod$TAC (written into om$catch)
		##    mod$TAC[,y] is total TAC, which is split into fishery and regions using om$catchsplits
		## 2. om$catch is split into sex proportions and stored in mod$real_catch
		## 3. Implementation error can be added (one form is added below for TAC)
		## Here: only static fleet dynamic (by using om$catchsplits to distribute TAC amongst gear types/fishery)
		## Possible additions for TAC: - Catch restrictions (eg changes in catch from year to year smaller than x ...)
		## 					  		   - Dynamic distribution of catch (change om$catchsplits for this)
		## If there is a TAC, use this:
		if(is.na(mod$TAC[,y]) == FALSE) { ##** better to if(!is.na(mod$TAC[,y]))
			om$catch[y,,,] <- as.numeric(mod$TAC[,y]) * om$catchsplits[y,,,]
			om$catch[y,,,] <- om$catch[y,,,] * exp(rnorm(1,mean=0,sd=om$TAC_implement_error))	# Implementation error
			om$catch[y,,,] <- round(om$catch[y,,,],0)
		}
		#### Seasonal Loop:
		for (ss in 1:om$n_seasons) {
			ssages <- om$ages + om$growth_props[ss]			# Ages for length and weight calculations
			# Move population numbers and tagging numbers from end of last season to start of this season (for season 1, see above)
			if (ss > 1) {
			  ## move the untagged population to the next season
				pop$n[,y,,ss,] <- pop$n[,y,,ss-1,]
				## move the tagged population
				if (sampling$pin_tagging) {
					if(sampling$pin_tag_Method == "Pool")
					  # Update season for pooled tags
						tag$tags[,y,,ss,,,] <- tag$tags[,y,,ss-1,,,]
					if(sampling$pin_tag_Method == "Ind") {
					  # Update season for individual tags
						liveT <- which(tag$iTags[,"Dead"] == 0)
						tag$iTags[liveT,"LastSeas"] <- tag$iTags[liveT,"LastSeas"] + 1
					}
				}
			}
			## 3. Maturation (if relevant for migration)
			##** placeholder
			## 4. Migration of pop numbers and tagged fish
			## implement different types of movement
			##* this looks like it overlaps the population array
			if(om$n_regions > 1) {
			  ## move the untagged population
			  switch(om$move_type,
			         array = {
			           stop("array based movement has not yet been implemented")
			           ##* this is array based movement
			           #pop$n[,y,,ss,] <- move_fish_array(n=pop$n[,y,,ss,,drop=FALSE],om=om)
			         },
			         dframe = {
			           ## this is vector based movement
			           pop$n[,y,,ss,] <- move_fish_dframe(n=pop$n[,y,,ss,,drop=FALSE],
			                                              move_rule=om$move_rules,
			                                              age_names=om$names_ages,
			                                              type="pop")
			         })
			  ## now move the tags
			  if(sampling$pin_tagging) {
			    if(sampling$pin_tag_Method == "Pool")
			      # Move pooled tags
			      tag$tags[,y,,ss,,,] <- move_fish_dframe(n=tag$tags[,y,,ss,,,,drop=FALSE],
			                                              move_rule=om$move_tag_rules,
			                                              age_names=om$names_ages,
			                                              type="tag")
			    ##* fork out this to another switch component
			    if(sampling$pin_tag_Method == "Ind" &  sum(tag$iTags[,"RelY"])) {
			      # Move individual tags
			      liveT <- which(tag$iTags[,"Dead"] == 0)
			      tag$iTags[liveT,] <- move_iTags(dat=tag$iTags[liveT,], om=om)
			    }
			  }
			}
			## 5. Natural Mortality: Remove half of natM
			half_natM <- exp(-0.5*pop$m[,y,,ss,])		# Natural Mortality
			pop$n[,y,,ss,] <- pop$n[,y,,ss,] * half_natM
			# Tagged fish: Remove half of natM (M of this year) of all (previously) tagged fish	with season-specific M
			if(sampling$pin_tagging){
				if(sampling$pin_tag_Method == "Pool"){ 	# NatM for pooled tags
					if(!sampling$pin_tag_Random) 		# Exact fractions
						tag$tags[,y,,ss,,,] <- sweep(tag$tags[,y,,ss,,,],c(1:5),half_natM,"*")
					if(sampling$pin_tag_Random) 		# Random binomial numbers due to small numbers
						tag$tags[,y,,ss,,,] <- tag_binom(dat1 = tag$tags[,y,,ss,,,], # Data
						                                 prob1 = half_natM, # Probabilities (age*sex*region)
						                                 dims = om$n_regions*om$n_years)	# Dimensions to expand half_natM
				}
				if(sampling$pin_tag_Method == "Ind" & sum(tag$iTags[,"RelY"]) > 0) {		# NatM for individual tags
					liveT <- which(tag$iTags[,"Dead"] == 0)
					tag$iTags[liveT,] <- iTags_M(dat = tag$iTags[liveT,],	# Data
					                             prob = half_natM)				# Probabilities (age*sex*region)
				}
			}
			#### For seasons with any catch:
			if(sum(om$catch[y,,ss,]) > 0) {
				## 6. Harvest rate
				for (ff in 1:length(om$fishery)) {
					for (rr in 1:om$n_regions) {		# Fishery-specific catch
						## For each region: Split total catch (om$catch) into sex proportions -> mod$real_catch
						if (length(om$sex) > 1) {
							TotBt  <- sum(fleet[[ff]]$landings_sel[,y,,ss,rr] * pop$n[,y,,ss,rr] * pop$wt[,y,,ss,rr])
							PartBt <- apply(fleet[[ff]]$landings_sel[,y,,ss,rr,drop=FALSE] * pop$n[,y,,ss,rr,drop=FALSE] * pop$wt[,y,,ss,rr,drop=FALSE],c(3),sum)
							mod$real_catch[[ff]][1,y,,ss,rr] <- om$catch[y,ff,ss,rr] * PartBt/TotBt
						} else {
							mod$real_catch[[ff]][1,y,1,ss,rr] <- om$catch[y,ff,ss,rr]
						}
						## Harvest rate by fishery (gear type) (= catch/ExplBt)
						mod$h_fishery[[ff]][1,y,,ss,rr] <- mod$real_catch[[ff]][,y,,ss,rr] /         		# Add 0.01 below to robustify
							(apply(fleet[[ff]]$landings_sel[,y,,ss,rr,drop=FALSE] * pop$n[,y,,ss,rr,drop=FALSE] * pop$wt[,y,,ss,rr,drop=FALSE],c(3),sum)+0.01)
						# Total H summed over fishery (gear types)
						mod$h_fishery_sum[1,y,,ss,rr] <- mod$h_fishery_sum[1,y,,ss,rr] + mod$h_fishery[[ff]][1,y,,ss,rr]
					}
				}
				## If total H > max H, reduce partial H (according to the weight of their contribution), no redistribution of catch
				## Possible change: Spatial fleet dynamics, eg higher catches where higher biomass or catch rates
				for (i in 1:length(om$sex)) {
					for (rr in 1:om$n_regions) {
						if (mod$h_fishery_sum[1,y,i,ss,rr] > mod$h_max) {
							h_fishery_diff <- mod$h_fishery_sum[,y,i,ss,rr] - mod$h_max	# Total difference
							for (ff in 1:length(om$fishery)) {    						# Weighted difference for each fishery
								weighted.diff <- (mod$h_fishery[[ff]][1,y,i,ss,rr]/mod$h_fishery_sum[,y,i,ss,rr]) * h_fishery_diff
								mod$h_fishery[[ff]][1,y,i,ss,rr] <- mod$h_fishery[[ff]][1,y,i,ss,rr] - weighted.diff
							}
							mod$h_fishery_sum[1,y,i,ss,rr] <- mod$h_max
						}
					}
				}
				## Calculate total H by age, landings and samples
				for (ff in 1:length(om$fishery)) {
					#ff <- 1
					for (rr in 1:om$n_regions) {
						## Total Harvest rate H by age: Sum of fishery-specific H * Sel by age
						## Used later to remove catch from pop$n
						mod$h_species[,y,,ss,rr] <- mod$h_species[,y,,ss,rr] +
										fleet[[ff]]$landings_sel[,y,,ss,rr] * mod$h_fishery[[ff]][1,y,,ss,rr]
						## Calculate catch-at-age by fishery (do not use h_species)
						fleet[[ff]]$landings_n[,y,,ss,rr] <- sweep(pop$n[,y,,ss,rr] * fleet[[ff]]$landings_sel[,y,,ss,rr],2,
						                                           mod$h_fishery[[ff]][1,y,,ss,rr],"*")
						## Calculate catch-at-length by fishery and summed over fisheries
						for (i in 1:length(om$sex)) {
							mod$landings_n_len[[ff]][,y,i,ss,rr] <- sample_lengths(ssages,sampling$len_classes,
							                                                      fleet[[ff]]$landings_n[,y,i,ss,rr],
							                                                      om$growth[[om$sex[i]]])
							mod$landings_n_len_sum[,y,i,ss,rr] <- mod$landings_n_len_sum[,y,i,ss,rr] + mod$landings_n_len[[ff]][,y,i,ss,rr]
						}
						## Actual catch taken using catch-by-age
						if (length(om$sex)  > 1){
							fleet[[ff]]$landings[,y,,ss,rr] <- apply(fleet[[ff]]$landings_n[,y,,ss,rr,drop=FALSE] * pop$wt[,y,,ss,rr,drop=FALSE],c(3),sum)
						}
						if (length(om$sex) == 1){
							fleet[[ff]]$landings[,y,,ss,rr] <- sum(fleet[[ff]]$landings_n[,y,,ss,rr] * pop$wt[,y,,ss,rr])
						}
						# Put back into mod$real_catch (om$catch is intended catch, mod$real_catch is actual catch)
						mod$real_catch[[ff]][,y,,ss,rr] <- fleet[[ff]]$landings[,y,,ss,rr]
						## Add observation error for catch:
						mod$obs_catch[[ff]][,y,,ss,rr] <- round(mod$real_catch[[ff]][,y,,ss,rr] * exp(rnorm(1,mean=0,sd=om$catch_obs_error)),0)
						## 7. Sampling
						## Sampling: age and length samples, absolute index of abundance (scaled by sampling$survey_q), cpue
						if (sampling$pin_sampling & sum(mod$real_catch[[ff]][,y,,ss,rr]) > 0) {
							obs[[ff]]$survey_index[,y,,ss,rr]	<- sample_survey(type="index", om=om, pop=pop, sampling=sampling, fleet=fleet, obs=obs, ff=ff,y=y,ss=ss,rr=rr)
							obs[[ff]]$survey_age_n[,y,,ss,rr] <- sample_survey(type="age", om=om, pop=pop, sampling=sampling, fleet=fleet, obs=obs, ff=ff,y=y,ss=ss,rr=rr)
							obs[[ff]]$survey_len_n[,y,,ss,rr] <- sample_survey(type="len", om=om, pop=pop, sampling=sampling, fleet=fleet, obs=obs, ff=ff,y=y,ss=ss,rr=rr)
							obs[[ff]]$age_sample_n[,y,,ss,rr] <- sample_catch(type="age", om=om, sampling=sampling, fleet=fleet, obs=obs, ff=ff,y=y,ss=ss,rr=rr)
							obs[[ff]]$len_sample_n[,y,,ss,rr] <- sample_catch(type="len", om=om, sampling=sampling, fleet=fleet, obs=obs, ff=ff,y=y,ss=ss,rr=rr)
							obs[[ff]]$cpue[,y,,ss,rr] <- sample_cpue(pop=pop, sampling=sampling, fleet=fleet, ff=ff, y=y, ss=ss, rr=rr)	# Relative index of abundance
						}
						## 8. Tagging and recapture
						if (sampling$pin_tagging & sum(mod$real_catch[[ff]][,y,,ss,rr]) > 0) {
							## Tagging
							#  if(sampling$pin_tag_N == "Fixed")  # Fixed number of tags released, no change needed
						  ##* extra brackets
							if (sampling$pin_tag_N == "Catch"){ 	# Tagging rate depending on catch and tagging rate
								obs[[ff]]$tag_N[,y,,ss,rr]  <- round(fleet[[ff]]$landings[,y,,ss,rr] * sampling$tag_rate[ff],0)
							}
							n.sel <- pop$n[,y,,ss,rr] * obs[[ff]]$tag_sel[,y,,ss,rr]	# N-at-age available for tagging by sex
							n.prop <- apply(n.sel,c(2),sum)/sum(n.sel)								# Proportion of N available by sex
							# tag_N is total number of tagged fish, not tagged fish by sex (as this depends on sex ratio)
							for (i in 1:length(om$sex)) {
								# Tag fish
								if (!sampling$pin_tag_Random) 		## Proportional: Multiply prop of n.sel by age with tag_N for each sex
									newtags <- (n.sel[,i]/sum(n.sel[,i])) * (obs[[ff]]$tag_N[,y,i,ss,rr]*n.prop[i])
								if (sampling$pin_tag_Random) 		## Multinomial random selection of tagged fish
									newtags <- rmultinom(n=1, size=obs[[ff]]$tag_N[,y,i,ss,rr]*n.prop[i],prob = n.sel[,i]/sum(n.sel[,i]))
								# Store released numbers-at-age of tagged fish in obs
								obs[[ff]]$tag_age_n[,y,i,ss,rr] <- newtags
								# Store released numbers-at-length of tagged fish in obs
								obs[[ff]]$tag_len_n[,y,i,ss,rr] <- sample_lengths(ssages,sampling$len_classes,obs[[ff]]$tag_age_n[,y,i,ss,rr],om$growth[[om$sex[i]]])
								# Account for release M & tag loss, and store tagged fish (pooled or individual tags)
								if(sampling$pin_tag_Method == "Pool") { 	# Pooled tags
									if (!sampling$pin_tag_Random) 		## Proportional tag-release mortality & tag loss
										newtags <- newtags*(1-sampling$tag_mort)
									if (sampling$pin_tag_Random) 		## Multinomial random tag-release mortality & tag loss
										newtags <- rbinom(n=om$n_ages, size=newtags, prob=(1-sampling$tag_mort))
									# Store pooled tags (sum over fishery)
									tag$tags[,y,i,ss,rr,y,rr] <- tag$tags[,y,i,ss,rr,y,rr] + newtags
								}
								if(sampling$pin_tag_Method == "Ind") { 	# Individual tags
									tag$iTags	<- iTags_release(dat=tag$iTags, newtags=newtags, y=y, i=i, ss=ss, rr=rr, om=om, sampling=sampling)
								}
							}
							## Recaptures
							## Recaptures of fish tagged in all previous years (current year is excluded):
							# Uses numbers-at-age for catch and population (before catch is being removed)
							prob.n <- fleet[[ff]]$landings_n[,y,,ss,rr]/pop$n[,y,,ss,rr]			# Prob of recapture by age class and sex
							prob.n[is.na(prob.n)]	<- 0
							if(sampling$pin_tag_Method == "Pool") { 	# Pooled tags
								if (!sampling$pin_tag_Random) 		# Proportional selection
									recap <- sweep(tag$tags[,y,,ss,rr,,],c(1,3),prob.n ,"*")
								if (sampling$pin_tag_Random) 		# Binomial random selection of recaptures - could also be multinomial?
									recap  <- tag_binom(dat1 = tag$tags[,y,,ss,rr,,],	# Data
									                    prob1 = prob.n,					# Probabilities
									                    dims = om$n_regions*om$n_years)
								if(om$n_regions==1) { recap[,,y] <- 0 } else { recap[,,y,] <- 0	}		# No recaptures of fish that were tagged in the same year
								# For each tagging year, sum over fishery
								tag$recaps[,y,,ss,rr,,] <- tag$recaps[,y,,ss,rr,,] + recap
								# Remove recaptured tags from fish.tagged
								tag$tags[,y,,ss,rr,,] <- tag$tags[,y,,ss,rr,,] - recap
								tag$tags[tag$tags[] < 0] <- 0	# Avoid negative numbers
							}
							if(sampling$pin_tag_Method == "Ind") { 	# Individual tags
									tag$iTags	<- iTags_recapture(dat=tag$iTags, prob=prob.n, y=y, ss=ss, rr=rr, om=om)
							}
						}	# Tag recapture loop
					}	# Region loop
				}	# Fishery loop
				## Calculate length from age for recaptures
				for (rr in 1:om$n_regions) {
					for (i in 1:length(om$sex)) {
						for (ytag in 1:y) {
							for (rtag in 1:om$n_regions) {
								# Recaptures: All recaptured in season ss and region rr of this year y from tagging events in past years ytag and regions rtag
								if (sum(tag$recaps[,y,i,ss,rr,ytag,rtag]) > 0) {
									tag$recaps_len[,y,i,ss,rr,ytag,rtag] <- sample_lengths(ssages,sampling$len_classes,tag$recaps[,y,i,ss,rr,ytag,rtag],om$growth[[om$sex[i]]])
									# Fish recaptured in year y and season ss
								}
							}
						}
					}
				}
			}	# If-loop for catch > 0
			## Store SSB and TotB (mid-season/year)
			for (rr in 1:om$n_regions) {
			  for (i in 1:length(om$sex)){
					mod$totB[,y,i,ss,rr] <- sum(pop$n[,y,i,ss,rr] * (1-0.5*mod$h_species[,y,i,ss,rr]) * pop$wt[,y,i,ss,rr])
  				mod$ssb[,y,i,ss,rr] <- sum(pop$n[,y,i,ss,rr] * (1-0.5*mod$h_species[,y,i,ss,rr]) * pop$wt[,y,i,ss,rr] * pop$fec[,y,i,ss,rr])
			  }
			}
			## Apply fishing mortality & second half of natM to fish stock: pop$n is now end-of-season numbers at age
			pop$n[,y,,ss,] <- pop$n[,y,,ss,] * (1-mod$h_species[,y,,ss,]) * exp(-0.5*pop$m[,y,,ss,])
			pop$n[,y,,ss,][pop$n[,y,,ss,] < 0] <- 0		# Safeguard that pop.n > 0
			# Tagged fish: Remove second half of natM (M of this year) of all (previously) tagged fish, and apply tag shedding rate (seasonal prop)
			half_natM_shedding <- exp(-0.5 * pop$m[,y,,ss,] - sampling$tag_shedding/om$n_seasons)
			if(sampling$pin_tag_Method == "Pool") { 	# Pooled tags
				if (!sampling$pin_tag_Random) 		# Exact fractions
					tag$tags[,y,,ss,,,] <- sweep(tag$tags[,y,,ss,,,], c(1:5), half_natM_shedding, "*")
				if (sampling$pin_tag_Random) 	# Random binomial numbers
					tag$tags[,y,,ss,,,] <- tag_binom(dat1  = tag$tags[,y,,ss,,,], # Data
													                 prob1 = half_natM_shedding)	# Probabilities
				tag$tags[,y,,ss,,,][tag$tags[,y,,ss,,,] < 0] <- 0	# Safeguard that tag$tags > 0
			}
			if(sampling$pin_tag_Method == "Ind" & sum(tag$iTags[,"RelY"]) > 0) { 	# Individual tags (with records)
				liveT <- which(tag$iTags[,"Dead"] == 0)
				tag$iTags[liveT,] <- 	iTags_M(dat = tag$iTags[liveT,], # Data
												              prob = half_natM_shedding)						# Probabilities (age*sex*region)
			}
		}	# Season loop
		#### Run CASAL assessment loop and TAC finder
		if(ctrl$pin_casal_assess == 1) {
			if(om$years[y] %in% ctrl$Assyr) {
				## Take catch of current year y as initial value for TAC search
				tac <- 0
				for (ff in 1:length(om$fishery))
					tac <- tac + round(as.numeric(sum(mod$real_catch[[ff]][,y,,,])), 0)
				## Run CASAL assessment (with adjustment for process error, with/without TAC finder)
				print("")
				print(paste("CASAL Assessment in year ", om$years[y], sep=""))
				## Update assessment parameters and add input data
				datass <- get_casal_data(Yr_current=om$years[y], datass=datass, om=om, sampling=sampling, obs=obs, tag=tag, mod=mod)
				# Update casal data if needed for running scenarios
				if(ctrl$pin_update_casal_data == 1) datass <- update_casal_data(datass=datass, Yr_current=om$years[y])
				## Create CASAL input files
				create_casal_file_est(params=datass, casal_path=ctrl$casal_path, skel_csl=ctrl$est_skel_csl, csl=ctrl$est_csl)
				create_casal_file_pop(params=datass, casal_path=ctrl$casal_path, skel_csl=ctrl$pop_skel_csl, csl=ctrl$pop_csl)
				create_casal_file_out(params=datass, casal_path=ctrl$casal_path, skel_csl=ctrl$output_skel_csl, csl=ctrl$output_csl)
				## Run assessment and TAC finder
			  mod$TAC[,y+1] <- run_complete_casal_assessment(ctrl=ctrl, om=om, datass=para[["ass"]], Yr_current=om$years[y],
			                                                 TAC_init=tac, intern = intern)
				## Store year-specific output.log and mdp.dat files
				nname <- unlist(strsplit(ctrl$output_log,"\\."))
				fname <- paste(ctrl$casal_path,nname[1],om$years[y],".",nname[2],sep="")
				file.copy(from=paste(ctrl$casal_path,ctrl$output_log,sep=""), to=fname, overwrite = TRUE)
				nname <- unlist(strsplit(ctrl$mpd_dat,"\\."))
				fname <- paste(ctrl$casal_path,nname[1],om$years[y],".",nname[2],sep="")
				file.copy(from=paste(ctrl$casal_path,ctrl$mpd_dat,sep=""), to=fname, overwrite = TRUE)
				nname <- unlist(strsplit(ctrl$proj_dat,"\\."))
				fname <- paste(ctrl$casal_path,nname[1],om$years[y],".",nname[2],sep="")
				file.copy(from=paste(ctrl$casal_path,ctrl$proj_dat,sep=""), to=fname, overwrite = TRUE)
			} else {
				## With previous assessment, but not an assessment year: roll over TAC
				if(om$years[y] > min(ctrl$Assyr)) mod$TAC[,y+1] <- mod$TAC[,y]
			}
		}
	}     # Annual loop # shouldn't this be end of annual loop?
	##* print statements slow things down
	print(paste("Total time for OM: ", round(as.vector(proc.time()-ptm)[3],1),sep=""))	# For time measurement
	## Return results to a list object
	obj <- list("pop" = pop,
	            "fleet" = fleet,
	            "obs" = obs,
	            "mod" = mod,
	            "tag" = tag)
	## return the results
	obj
}
