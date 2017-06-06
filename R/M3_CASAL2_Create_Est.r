#####	Create estimation.csl file for CASAL Assessment

#' Create CASAL estimation.csl file
#'
#' Create CASAL estimation.csl file
#' @inheritParams casal
#' @export
create_casal_file_est <- function(params, casal_path, skel_csl, csl) {
	#### This funtion sets up a casal_population file (perhaps estimation file) for an assessment
	# params: Object with input parameters for CASAL assessment
	# params    	<- para[["ass"]]
	# casal_path 	<- ""
	# skel_csl  	<- "casal_estimation_skel.csl"
	# csl 			<- "casal_estimation.csl"
	# Extract casal parameters from skeleton file
	casalest <- casal::extract.csl.file(paste(casal_path,skel_csl,sep=""))
	## Get assessment year
	ycurr <- which(params$years==params$Yr_current)
	#### 0. Prepare csl file (for sequence, see Manual Chapter 14.1.2 The input parameter files)
	# Remove some list elements (that relate to individual fisheries or tagging data)
	casalest[["profile[1]"]] <- NULL
	casalest[["relative_abundance[f1_cpue]"]] <- NULL
	casalest[["relative_numbers_at[Survey1_Age]"]] <- NULL
	casalest[["relative_numbers_at[Survey2_Size]"]] <- NULL
	casalest[["catch_at[Catch_f1Age]"]] <- NULL
	casalest[["catch_at[Catch_f1Size]"]] <- NULL
	casalest[["estimate[1]"]] <- NULL
	casalest[["catch_limit_penalty[1]"]] <- NULL
	casalest[["vector_average_penalty[1]"]] <- NULL
	#### 1. Estimation
	casalest[["estimator"]]$value <- params$estimator
	casalest[["max_iters"]]$value <- params$max_iters
	casalest[["max_evals"]]$value <- params$max_evals
	casalest[["grad_tol"]]$value <- params$grad_tol
	casalest[["MCMC"]]$command <- "MCMC"
	casalest[["MCMC"]]$value <- character(0)
	casalest[["MCMC"]]$start <- params$MCMC_start
	casalest[["MCMC"]]$length <- params$MCMC_length
	casalest[["MCMC"]]$keep <- params$MCMC_keep
	casalest[["MCMC"]]$stepsize <- params$MCMC_stepsize
	casalest[["MCMC"]]$adaptive_stepsize <- params$MCMC_adaptive_stepsize
	casalest[["MCMC"]]$adapt_at <- params$MCMC_adapt_at
	casalest[["MCMC"]]$burn_in <- params$MCMC_burn_in
	#### 2. Estimation profiling
	casalest[["profile[1]"]]$command <- "profile"
	casalest[["profile[1]"]]$value <- character(0)
	casalest[["profile[1]"]]$parameter <- params$profile_parameter
	casalest[["profile[1]"]]$n <- params$profile_n
	casalest[["profile[1]"]]$l <- params$profile_l
	casalest[["profile[1]"]]$u <- params$profile_u
	#### 3. Relative index of abundance: CPUE
	list_abund <- NULL
	for (ff in 1:length(params$list_fishery)) {
		fishery <- paste(params$list_fishery[ff],"_cpue",sep="")
		fisheryN <- paste("relative_abundance[",fishery,"]",sep="")
		fish <- params$Fish[ff,]
		if (sum(params$catch[[params$list_fishery[ff]]]) > 0 & sum(params$sample_years["cpue",1:ycurr,ff]) > 0) {	# Catch > 0 and Sampling = 1
			# CPUE Index for years with sampling = 1 (obsyy) and catch > 0 (catchyy)
			yy1 <- params$sample_years["cpue",1:ycurr,fish["Fishery"]]
			obsyy <- dimnames(params$sample_years[,yy1 > 0,])$year
			catchbyyear <- apply(params$catch[[params$list_fishery[ff]]],2,sum)
			catchyy <- names(catchbyyear[catchbyyear > 0])
			obsyy <- obsyy[obsyy %in% catchyy]
			dat <- apply(params$cpue_index[[fish["Fishery"]]][,obsyy,,,],c(1),sum)	# Sum cpue_index (is catch-weighted by sex)
			dat[dat==0] <- 0.01											# Robustify against lognormal distribution
			datcv <- params$cpue_cv[[fish["Fishery"]]][,obsyy,1,,]	# Take cv from first flunit
			# Data for estimation.csl file
			casalest[[fisheryN]]$command <- "relative_abundance"
			casalest[[fisheryN]]$value <- fishery
			casalest[[fisheryN]]$biomass <- "True"
			casalest[[fisheryN]]$area <- as.vector(fish["Region"])		# Region
			casalest[[fisheryN]]$step <- as.vector(fish["Season"])        # Season
			casalest[[fisheryN]]$q <- as.vector(fish["q"])     		# Catchability
			casalest[[fisheryN]]$ogive <- as.vector(fish["Sel"])  			# Selectivity
			casalest[[fisheryN]]$proportion_mortality <- params$proportion_mortality	# Prop mortality before observation
			casalest[[fisheryN]]$years <- c(names(dat))     				# List of years with observations
			for (yy in 1:length(obsyy)) casalest[[fisheryN]][names(dat)[yy]] <- round(dat[yy],4)
			for (yy in 1:length(obsyy)) casalest[[fisheryN]][paste("cv_",names(dat)[yy],sep="")] <- round(datcv[yy],4)
			casalest[[fisheryN]]$dist <- params$cpue_cv_dist
			casalest[[fisheryN]]$cv_process_error	<- params$cpue_cv_process_error
			# List used below in @estimate
			list_abund <- c(list_abund,fishery)
		}
	}
	#### 4. Survey index
	list_abund <- NULL
	for (ff in 1:length(params$list_fishery)) {
		fishery <- paste(params$list_fishery[ff],"_survIndex",sep="")
		fisheryN <- paste("relative_abundance[",fishery,"]",sep="")
		fish <- params$Fish[ff,]
		if (sum(params$catch[[params$list_fishery[ff]]]) > 0 & sum(params$sample_years["survey",1:ycurr,ff]) > 0) {		# Catch > 0 and Sampling = 1
			# Abundance Index for years sampled obsyy (summed over sex, cv taken from sex[1])
			yy1 <- params$sample_years["survey",1:ycurr,fish["Fishery"]]
			obsyy <- dimnames(params$sample_years[,yy1 > 0,])$year
			catchbyyear <- apply(params$catch[[params$list_fishery[ff]]],2,sum)
			catchyy <- names(catchbyyear[catchbyyear > 0])
			obsyy <- obsyy[obsyy %in% catchyy]
			dat <- apply(params$survey_index[[fish["Fishery"]]][,obsyy,,,],c(1),sum)
			dat[dat==0] <- 0.01													# Robustify against lognormal distribution
			datcv <- params$survey_index_var[[fish["Fishery"]]][,obsyy,1,,]	# Take cv from first flunit
			# Data for estimation.csl file (area, step, q and selectivity defined through fishery)
			casalest[[fisheryN]]$command <- "relative_abundance"
			casalest[[fisheryN]]$value <- fishery
			casalest[[fisheryN]]$biomass <- "True"
			casalest[[fisheryN]]$area <- as.vector(fish["Region"])		# Region
			casalest[[fisheryN]]$step <- as.vector(fish["Season"])        # Season
			casalest[[fisheryN]]$q <- as.vector(fish["q"])     		# Catchability
			casalest[[fisheryN]]$ogive <- as.vector(fish["Sel"])  			# Selectivity
			casalest[[fisheryN]]$proportion_mortality <- params$proportion_mortality	# Prop mortality before observation
			casalest[[fisheryN]]$years <- c(names(dat))    				# List of years with observations
			for (yy in 1:length(obsyy)) casalest[[fisheryN]][obsyy[yy]] <- round(dat[yy],1)
			for (yy in 1:length(obsyy)) casalest[[fisheryN]][paste("cv_",obsyy[yy],sep="")] <- datcv[yy]
			casalest[[fisheryN]]$dist <- params$survey_dist
			casalest[[fisheryN]]$cv_process_error	<- params$survey_cv_process_error
			# List used below in @estimate
			list_abund <- c(list_abund,fishery)
		}
	}
	#### 5. Survey Catch-at-age
	list_survA <- NULL
	for (ff in 1:length(params$list_fishery)) {
		fishery <- paste(params$list_fishery[ff],"_survA",sep="")
		fisheryN <- paste("relative_numbers_at[",fishery,"]",sep="")
		fish <- params$Fish[ff,]
		if (sum(params$catch[[params$list_fishery[ff]]]) > 0 & sum(params$sample_years["surveyage",1:ycurr,ff]) > 0) {		# Catch > 0 and Sampling = 1
			# Data (sum over sex if single-sex)
			yy1 <- params$sample_years["surveyage",1:ycurr,fish["Fishery"]]
			obsyy <- dimnames(params$sample_years[,yy1 > 0,])$year
			catchbyyear <- apply(params$catch[[params$list_fishery[ff]]],2,sum)
			catchyy <- names(catchbyyear[catchbyyear > 0])
			obsyy <- obsyy[obsyy %in% catchyy]
			dat <- params$survey_age_n[[fish["Fishery"]]][,obsyy,,,]
			dat[is.na(dat)] <- 0
			trun <- which(dimnames(dat)$age==max(as.integer(dimnames(dat)$age[apply(dat,c(1),sum)>0]))) # Range with obs > 0 (summed over all years)
			dat <- dat[1:trun,,]
			datcv <- params$survey_age_var[[fish["Fishery"]]][1:trun,obsyy,,,]
			if (params$by_sex == 0) {		# Single sex
				dat <- apply(dat,c(1,2),sum)
				datcv <- apply(datcv,c(1,2),mean)
			}
			dat[dat==0] <- 1		# Robustify against lognormal LL
			# Store in casalest
			casalest[[fisheryN]]$command <- "relative_numbers_at"
			casalest[[fisheryN]]$value <- fishery
			casalest[[fisheryN]]$area <- as.vector(fish["Region"])		# Region
			casalest[[fisheryN]]$step <- as.vector(fish["Season"])        # Season
			casalest[[fisheryN]]$q <- as.vector(fish["q"])  	     	# Catchability
			casalest[[fisheryN]]$ogive <- as.vector(fish["Sel"])			# Selectivity
			casalest[[fisheryN]]$proportion_mortality <- params$proportion_mortality	# Prop mortality before observation
			casalest[[fisheryN]]$at_size <- "False"
			casalest[[fisheryN]]$sexed <- if(params$by_sex == 0) "False" else "True"
			casalest[[fisheryN]]$plus_group <- params$age_plus_group
			casalest[[fisheryN]]$min_class <- params$ages[1]
			casalest[[fisheryN]]$max_class <- params$ages[trun]
			casalest[[fisheryN]]$years <- c(dimnames(dat)$year)        	# List of years with observations
			if (params$by_sex == 0) {
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <- paste(dat[,dimnames(dat)$year==obsyy[yy]],collapse=" ")
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][paste("cvs_",obsyy[yy],sep="")] <- paste(datcv[,dimnames(datcv)$year==obsyy[yy]],collapse=" ")
			} else {			# With sex partition: CASAL requires first males then females (in the same line)
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <-
						paste(paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="m"],collapse=" "),
							  paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="f"],collapse=" "))
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][paste("cvs_",obsyy[yy],sep="")] <-
						paste(paste(datcv[,dimnames(datcv)$year==obsyy[yy],dimnames(datcv)$sex=="m"],collapse=" "),
							  paste(datcv[,dimnames(datcv)$year==obsyy[yy],dimnames(datcv)$sex=="f"],collapse=" "))
			}
			casalest[[fisheryN]]$dist <- params$survey_dist
			casalest[[fisheryN]]$cv_process_error <- params$survey_cv_process_error
			# List used below in @estimate
			list_survA <- c(list_survA, fishery)
		}
	}
	#### 6. Survey Catch-at-length
	list_survL <- NULL
	for (ff in 1:length(params$list_fishery)) {
		fishery <- paste(params$list_fishery[ff],"_survL",sep="")
		fisheryN <- paste("relative_numbers_at[",fishery,"]",sep="")
		fish <- params$Fish[ff,]
		if (sum(params$catch[[params$list_fishery[ff]]]) > 0 & sum(params$sample_years["surveylen",1:ycurr,ff]) > 0) {		# Catch > 0 and Sampling = 1
			# Data (sum over sex if single-sex)
			yy1 <- params$sample_years["surveylen",1:ycurr,fish["Fishery"]]
			obsyy <- dimnames(params$sample_years[,yy1 > 0,])$year
			catchbyyear <- apply(params$catch[[params$list_fishery[ff]]],2,sum)
			catchyy <- names(catchbyyear[catchbyyear > 0])
			obsyy <- obsyy[obsyy %in% catchyy]
			dat <- params$survey_len_n[[fish["Fishery"]]][,obsyy,,,]
			dat[is.na(dat)] <- 0
			trun <- which(dimnames(dat)$len==max(as.integer(dimnames(dat)$len[apply(dat,c(1),sum)>0]))) # Range with obs > 0 (summed over all years)
			dat <- dat[1:trun,,]
			datcv <- params$survey_len_var[[fish["Fishery"]]][1:trun,obsyy,,,]
			if (params$by_sex == 0) {		# Single sex
				dat <- apply(dat,c(1,2),sum)
				datcv <- apply(datcv,c(1,2),mean)
			}
			dat[dat==0] <- 1		# Robustify against lognormal LL
			# Data for estimation.csl file
			casalest[[fisheryN]]$command <- "relative_numbers_at"
			casalest[[fisheryN]]$value <- fishery
			casalest[[fisheryN]]$area <- as.vector(fish["Region"])		# Region
			casalest[[fisheryN]]$step <- as.vector(fish["Season"])        # Season
			casalest[[fisheryN]]$q <- as.vector(fish["q"])  	     	# Catchability
			casalest[[fisheryN]]$ogive <- as.vector(fish["Sel"])			# Selectivity
			casalest[[fisheryN]]$proportion_mortality <- params$proportion_mortality	# Prop mortality before observation
			casalest[[fisheryN]]$at_size <- "True"
			casalest[[fisheryN]]$sexed <- if(params$by_sex == 0) "False" else "True"
			casalest[[fisheryN]]$plus_group <- params$len_plus_group
			casalest[[fisheryN]]$years <- c(dimnames(dat)$year) 		   	# List of years with observations
			casalest[[fisheryN]]$class_mins <- if(params$len_plus_group == "False") params$class_mins[1:(trun+1)] else params$class_mins[1:trun]
			if (params$by_sex == 0) {
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <- paste(dat[,dimnames(dat)$year==obsyy[yy]],collapse=" ")
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][paste("cvs_",obsyy[yy],sep="")] <- paste(datcv[,dimnames(datcv)$year==obsyy[yy]],collapse=" ")
			} else {			# With sex partition: CASAL requires first males then females (in the same line)
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <-
						paste(paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="m"],collapse=" "),
							  paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="f"],collapse=" "))
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][paste("cvs_",obsyy[yy],sep="")] <-
						paste(paste(datcv[,dimnames(datcv)$year==obsyy[yy],dimnames(datcv)$sex=="m"],collapse=" "),
							  paste(datcv[,dimnames(datcv)$year==obsyy[yy],dimnames(datcv)$sex=="f"],collapse=" "))
			}
			casalest[[fisheryN]]$dist <- params$survey_dist
			casalest[[fisheryN]]$cv_process_error <- params$survey_cv_process_error
			# List used below in @estimate
			list_survL <- c(list_survL,fishery)
		}
	}
	#### 7. Catch-at-age
	list_catchA <- NULL
	for (ff in 1:length(params$list_fishery)) {
		fishery <- paste(params$list_fishery[ff],"_catchA",sep="")
		fisheryN <- paste("catch_at[",fishery,"]",sep="")
		fish <- params$Fish[ff,]
		if (sum(params$catch[[params$list_fishery[ff]]]) > 0 & sum(params$sample_years["catchage",1:ycurr,ff]) > 0) {		# Catch > 0 and Sampling = 1
			# Data (sum over sex, datN is total numbers)
			yy1 <- params$sample_years["catchage",1:ycurr,fish["Fishery"]]
			obsyy <- dimnames(params$sample_years[,yy1 > 0,])$year
			catchbyyear <- apply(params$catch[[params$list_fishery[ff]]],2,sum)
			catchyy <- names(catchbyyear[catchbyyear > 0])
			obsyy <- obsyy[obsyy %in% catchyy]
			dat <- params$catch_age[[fish["Fishery"]]][,obsyy,,,]
			dat[is.na(dat)] <- 0
			trun <- which(dimnames(dat)$age==max(as.integer(dimnames(dat)$age[apply(dat,c(1),sum)>0]))) # Range with obs > 0 (summed over all years)
			dat	<- dat[1:trun,,]
			datN <- apply(dat,2,sum)
			if (params$by_sex == 0) dat	<- apply(dat,c(1,2),sum)				# Single-sex
			dat <- sweep(dat,2,datN,"/")	# Proportions
			dat[dat==Inf] <- 0						# Where sample size = 0, replace Inf with 0
			# dat[dat==0] <- 0.01						# Robustify against lognormal LL
			dat <- round(dat,4)
			# Data for estimation.csl file (area, step, q and selectivity defined through fishery)
			casalest[[fisheryN]]$command <- "catch_at"
			casalest[[fisheryN]]$value <- fishery						# fishery_catchA
			casalest[[fisheryN]]$fishery <- params$list_fishery[ff]		# fishery
			casalest[[fisheryN]]$at_size <- "False"
			casalest[[fisheryN]]$sexed <- if(params$by_sex == 0) "False" else "True"
			casalest[[fisheryN]]$sum_to_one <- "False"
			casalest[[fisheryN]]$plus_group <- params$age_plus_group
			casalest[[fisheryN]]$min_class <- params$ages[1]
			casalest[[fisheryN]]$max_class <- params$ages[trun]
			casalest[[fisheryN]]$years <- c(dimnames(dat)$year)     	# List of years with observations
			if(params$by_sex == 0) {
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <- paste(dat[,dimnames(dat)$year==obsyy[yy]],collapse=" ")
			} else {			# With sex partition: CASAL requires first males then females (in the same line)
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <-
						paste(paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="m"],collapse=" "),
							  paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="f"],collapse=" "))
			}
			for (yy in 1:length(obsyy)) {
				if(params$pin_diff_N_in_assessment == 1) {		# Adjust sample sizes if required
					casalest[[fisheryN]][paste("N_",obsyy[yy],sep="")] <- params$Fish[ff,"catchage_N"] # params$catchage_N[[ff]]
				} else {										# Otherwise take sample sizes from original samples in OM
					casalest[[fisheryN]][paste("N_",obsyy[yy],sep="")] <- datN[yy]
				}
			}
			casalest[[fisheryN]]$dist <- params$catch_at_dist
			casalest[[fisheryN]]$r <- params$catch_at_r
			# List used below in @estimate
			list_catchA <- c(list_catchA,fishery)
		}
	}
	#### 8. Catch-at-length
	list_catchL <- NULL
	for (ff in 1:length(params$list_fishery)) {
		fishery <- paste(params$list_fishery[ff],"_catchL",sep="")
		fisheryN <- paste("catch_at[",fishery,"]",sep="")
		fish <- params$Fish[ff,]
		if (sum(params$catch[[params$list_fishery[ff]]]) > 0 & sum(params$sample_years["catchlen",1:ycurr,ff]) > 0) {		# Catch > 0 and Sampling = 1
			# Data (sum over sex, N is total numbers)
			yy1 <- params$sample_years["catchage",1:ycurr,fish["Fishery"]]
			obsyy <- dimnames(params$sample_years[,yy1 > 0,])$year
			catchbyyear <- apply(params$catch[[params$list_fishery[ff]]],2,sum)
			catchyy <- names(catchbyyear[catchbyyear > 0])
			obsyy <- obsyy[obsyy %in% catchyy]
			dat <- params$catch_len[[fish["Fishery"]]][,obsyy,,,]
			dat[is.na(dat)] <- 0
			trun <- which(dimnames(dat)$len==max(as.integer(dimnames(dat)$len[apply(dat,c(1),sum)>0]))) # Range with obs > 0 (summed over all years)
			dat <- dat[1:trun,,]
			datN <- apply(dat,2,sum)
			if (params$by_sex == 0) dat	<- apply(dat,c(1,2),sum)				# Single-sex
			dat <- sweep(dat,2,datN,"/")	# Proportions
			dat[dat==Inf] <- 0						# Where sample size = 0, replace Inf with 0
			# dat[dat==0] <- 0.01						# Robustify against lognormal LL
			dat <- round(dat,4)
			# Data for estimation.csl file (no area, step, q or selectivity (defined through fishery)
			casalest[[fisheryN]]$command <- "catch_at"
			casalest[[fisheryN]]$value <- fishery						# fishery_catchL
			casalest[[fisheryN]]$fishery <- params$list_fishery[ff]		# fishery
			casalest[[fisheryN]]$at_size <- "True"
			casalest[[fisheryN]]$sexed <- if(params$by_sex == 0) "False" else "True"
			casalest[[fisheryN]]$sum_to_one <- "False"
			casalest[[fisheryN]]$plus_group <- params$len_plus_group
			casalest[[fisheryN]]$years <- c(dimnames(dat)$year)       	# List of years with observations
			casalest[[fisheryN]]$class_mins <- if(params$len_plus_group == "False") params$class_mins[1:(trun+1)] else params$class_mins[1:trun]
			if(params$by_sex == 0) {
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <- paste(dat[,dimnames(dat)$year==obsyy[yy]],collapse=" ")
			} else {			# With sex partition: CASAL requires first males then females (in the same line)
				for (yy in 1:length(obsyy))
					casalest[[fisheryN]][obsyy[yy]] <-
						paste(paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="m"],collapse=" "),
							  paste(dat[,dimnames(dat)$year==obsyy[yy],dimnames(dat)$sex=="f"],collapse=" "))
			}
			for (yy in 1:length(obsyy)) {
				if(params$pin_diff_N_in_assessment == 1) {		# Adjust sample sizes if required
					casalest[[fisheryN]][paste("N_",obsyy[yy],sep="")] <- params$Fish[ff,"catchlen_N"]# params$catchage_N[[ff]]
				} else {										# Otherwise take sample sizes from original samples in OM
					casalest[[fisheryN]][paste("N_",obsyy[yy],sep="")] <- datN[yy]
				}
			}
			casalest[[fisheryN]]$dist <- params$catch_at_dist
			casalest[[fisheryN]]$r <- params$catch_at_r
			# List used below in @estimate
			list_catchL <- c(list_catchL,fishery)
		}
	}
	#### 9. Tag-recaptures
	## Recaptures can be season and fishery-specific (with ogive = fishery selectivity)
	##     or pooled across seasons and fisheries (without ogive/selectivity: default=no selectivity)
	## To limit the number of tag partitions, recaptures are pooled for all fisheries and fishing seasons
	## If changes are needed, change OM and add loops for params$list_fishery and seasons
	list_tagrecap <- NULL
	ttag <- params$sample_years["tagging",1:ycurr,,drop=FALSE]
	if(sum(ttag)>0) {
		## tag_names: Vector of tagging release partitions, links releases with recaptures (combination of year * region)
		revents <- params$Rec_events
		## loop over the release events
		for (ff in 1:nrow(revents)) {
			tagN <- paste("tag_recapture[",revents[ff,"RecName"],"]",sep="")
			dat <- params$Rec_data[[revents[ff,"RecName"]]]
			scannedN <- params$scannedN[[revents[ff,"RecName"]]]
			ryear <- dimnames(dat)$year
			if (sum(dat) > 0) {
				# Write casalest: Only tag_sampling_type == "size" implemented so far
				casalest[[tagN]]$command <- "tag_recapture"
				casalest[[tagN]]$value <- revents[ff,"RecName"]
				casalest[[tagN]]$tag_name <- revents[ff,"TagName"]
				casalest[[tagN]]$area <- revents[ff,"RecRegion"]
				casalest[[tagN]]$step <- params$Fish[1,"Season"]         	# Season: Take Season from first fishery
				casalest[[tagN]]$proportion_mortality <- params$proportion_mortality		# Prop mortality before observation
				casalest[[tagN]]$detection_probability <- params$tag_detection_probability
				casalest[[tagN]]$sample <- params$tag_sampling_type
				casalest[[tagN]]$years <- ryear
				casalest[[tagN]]$plus_group <- params$len_plus_group
				casalest[[tagN]]$class_mins <- params$class_mins
				if (params$by_sex == 0) {
					for (y in 1:length(ryear))
						casalest[[tagN]][paste("recaptured_",ryear[y],sep="")] <- paste(dat[,y],collapse=" ")
					for (y in 1:length(ryear))
						casalest[[tagN]][paste("scanned_",ryear[y],sep="")] <- paste(scannedN[,y],collapse=" ")
				} else {	# With sex partition
					for (y in 1:length(ryear)) {
						casalest[[tagN]][paste("recaptured_female_",ryear[y],sep="")]	<- paste(dat[,y,dimnames(dat)$sex=="f"],collapse=" ")
						casalest[[tagN]][paste("recaptured_male_",ryear[y],sep="")] <- paste(dat[,y,dimnames(dat)$sex=="m"],collapse=" ")
					}
					for (y in 1:length(ryear)) {
						casalest[[tagN]][paste("scanned_female_",ryear[y],sep="")] <- paste(scannedN[,y,dimnames(scannedN)$sex=="f"],collapse=" ")
						casalest[[tagN]][paste("scanned_male_",ryear[y],sep="")] <- paste(scannedN[,y,dimnames(scannedN)$sex=="m"],collapse=" ")
					}
				}
				casalest[[tagN]]$do_bootstrap <- params$tag_do_bootstrap
				casalest[[tagN]]$r <- params$tag_r
				casalest[[tagN]]$dispersion <- params$tag_dispersion
				# List used below in @estimate
				list_tagrecap <- c(list_tagrecap,fishery)
			}
		}
	}
	#### Relativity constants (Estimation for q)
	# Ageing error
	casalest[["ageing_error"]]$type <- params$ageing_error[1]
	casalest[["ageing_error"]]$c <- params$ageing_error[2]
	# Count of '@estimate' list elements
	estN <- 0
	# Catchabilities
	casalest[["q_method"]]$value <- params$q_method		# @q_method nuisance
	for (ee in 1:length(params$qq_names)) {					#  qq_names hold unique catchability names (<> list_q)
		estN <- estN + 1
		estimateN <- paste("estimate[",estN,"]",sep="")
		casalest[[estimateN]]$command <- "estimate"
		casalest[[estimateN]]$value <- character(0)
		casalest[[estimateN]]$parameter <- paste("q[",params$qq_names[ee],"].q",sep="")
		casalest[[estimateN]]$lower_bound <- params$qqvalues[[ee]][1]
		casalest[[estimateN]]$upper_bound <- params$qqvalues[[ee]][2]
		casalest[[estimateN]]$prior <- params$qqvalues[[ee]][3]
		if (params$qqvalues[[ee]][3] == "lognormal") {
			casalest[[estimateN]]$mu <- params$qqvalues[[ee]][4]
			casalest[[estimateN]]$cv <- params$qqvalues[[ee]][5]
		}
		casalest[[estimateN]]$phase <- 1
	}
	#### Estimation for free parameters
	for (ll in 1:4) {
		if (ll == 1) freeparam <- params$estim_initialization.B0
		if (ll == 2) freeparam <- params$estim_size_at_age.cv
		if (ll == 3) freeparam <- params$estim_natural_mortality.all
		if (ll == 4) freeparam <- params$estim_natural_mortality.ogive_all
		if(freeparam[1] == 1) {
			estN <- estN + 1
			estimateN <- paste("estimate[",estN,"]",sep="")
			casalest[[estimateN]]$command <- "estimate"
			casalest[[estimateN]]$value <- character(0)
			casalest[[estimateN]]$parameter <- freeparam[5]
			casalest[[estimateN]]$value <- character(0)
			casalest[[estimateN]]$lower_bound <- freeparam[2]
			casalest[[estimateN]]$upper_bound <- freeparam[3]
			casalest[[estimateN]]$prior <- freeparam[4]
			if(casalest[[estimateN]]$prior == "lognormal") {
				casalest[[estimateN]]$mu <- freeparam[6]
				casalest[[estimateN]]$cv <- freeparam[7]
			}
		}
	}
	if(params$estim_recruitment.YCS[1] == 1) {
		estN <- estN + 1
		estimateN <- paste("estimate[",estN,"]",sep="")
		casalest[[estimateN]]$command <- "estimate"
		casalest[[estimateN]]$value <- character(0)
		casalest[[estimateN]]$parameter <- params$estim_recruitment.YCS[[5]]
		casalest[[estimateN]]$value <- character(0)
		casalest[[estimateN]]$lower_bound <- params$estim_recruitment.YCS[[2]]
		casalest[[estimateN]]$upper_bound <- params$estim_recruitment.YCS[[3]]
		casalest[[estimateN]]$prior <- params$estim_recruitment.YCS[[4]]
		if(casalest[[estimateN]]$prior == "lognormal") {
			casalest[[estimateN]]$mu <- params$estim_recruitment.YCS[[6]]
			casalest[[estimateN]]$cv <- params$estim_recruitment.YCS[[7]]
		}
	}
	# cv_process_error
	for (ll in 1:3) {
		if (ll == 1)  {
			llist <- list_abund
			freeparam <- params$estim_abund_cv_process_error
			prefix <- "relative_abundance"
		}
		if (ll == 2) {
			llist <- list_survA
			freeparam <- params$estim_survA_cv_process_error
			prefix <- "relative_numbers_at"
		}
		if (ll == 3) {
			llist <- list_survL
			freeparam <- params$estim_survS_cv_process_error
			prefix <- "relative_numbers_at"
		}
		if(freeparam[1] == 1 & length(llist) > 0) {
			for (ee in 1: length(llist)){
				estN <- estN + 1
				estimateN <- paste("estimate[",estN,"]",sep="")
				casalest[[estimateN]]$command <- "estimate"
				casalest[[estimateN]]$value <- character(0)
				casalest[[estimateN]]$parameter <- paste(prefix,"[",llist[ee],"].cv_process_error",sep="")
				casalest[[estimateN]]$lower_bound <- freeparam[2]
				casalest[[estimateN]]$upper_bound <- freeparam[3]
				casalest[[estimateN]]$prior <- freeparam[4]
			}
		}
	}
	##* This creates the selectivity in the estimation file
	# Selectivities
	# estimate_selectivity holds selectivity names to estimate (<> list_sel)
	if(!is.null(params$estimate_selectivity)){
	  for (ee in 1:length(params$estimate_selectivity)) {
		  estN <- estN + 1
		  estimateN <- paste("estimate[",estN,"]",sep="")
		  casalest[[estimateN]]$command <- "estimate"
		  casalest[[estimateN]]$value <- character(0)
		  casalest[[estimateN]]$parameter <- paste("selectivity[",params$estimate_selectivity[ee],"].all",sep="")
		  if (params$by_sex==0) {
			  casalest[[estimateN]]$lower_bound <- params$est_selN_all[[ee]][[1]]
			  casalest[[estimateN]]$upper_bound <- params$est_selN_all[[ee]][[2]]
			  casalest[[estimateN]]$prior <- params$est_selN_all[[ee]][[3]]
		  }
	  }
	}
	## Migration
	##** add if(n_migrations>0)...
	if(params$migration_est_pin	== TRUE) {
		for (ee in 1:length(params$migration_names)) {
			estN <- estN + 1
			estimateN <- paste("estimate[",estN,"]",sep="")
			casalest[[estimateN]]$command <- "estimate"
			casalest[[estimateN]]$value	<- character(0)
			casalest[[estimateN]]$parameter <- paste("migration[",params$migration_names[ee],"].rates_all",sep="")
			casalest[[estimateN]]$lower_bound <- params$migration_rates_all_low
			casalest[[estimateN]]$upper_bound <- params$migration_rates_all_upp
			casalest[[estimateN]]$prior <- params$migration_rates_all_prior
		}
	}
	#### Penalties for all observations
	for (ee in 1: length(params$list_fishery)){
		estimateN <- paste("catch_limit_penalty[",ee,"]",sep="")
		casalest[[estimateN]]$command <- "catch_limit_penalty"
		casalest[[estimateN]]$value <- character(0)
		casalest[[estimateN]]$label <- paste(params$list_fishery[ee],"_Penalty",sep="")
		casalest[[estimateN]]$fishery <- params$list_fishery[ee]
		casalest[[estimateN]]$log_scale <- params$catch_limit_penalty_log_scale
		casalest[[estimateN]]$multiplier <- params$catch_limit_penalty_multiplier
	}
	# Tagging penalties for each tag event (corresponds to @tag command in the population.csl)
	if (sum(ttag) > 0) {
		count <- 0
		tag_names <- params$tag_names
		if(nrow(tag_names)>0) {
			## Same code as in population.csl
			for (ff in 1:nrow(tag_names)) {
				count <- count + 1
				tevent <- tag_names[ff,"TagEvent"]
				estimateN <- paste("fish_tagged_penalty[",count,"]",sep="")
				casalest[[estimateN]]$command <- "fish_tagged_penalty"
				casalest[[estimateN]]$value	<- character(0)
				casalest[[estimateN]]$label <- paste("Penalty_",tevent,sep="")
				casalest[[estimateN]]$tagging_episode <- tevent
				casalest[[estimateN]]$multiplier <- params$fish_tagged_penalty_multiplier
			}
		}
	}
	# Vector_average_penalty
	if(params$vector_average_penalty[1] == 1) {
		casalest[["vector_average_penalty"]]$command <- "vector_average_penalty"
		casalest[["vector_average_penalty"]]$value <- character(0)
		casalest[["vector_average_penalty"]]$label <- params$vector_average_penalty[2]
		casalest[["vector_average_penalty"]]$vector <- params$vector_average_penalty[3]
		casalest[["vector_average_penalty"]]$k <- params$vector_average_penalty[4]
		casalest[["vector_average_penalty"]]$multiplier <- params$vector_average_penalty[5]
	}
	# write the casal estimation file
	casal::write.csl.file(casalest, paste(casal_path, csl,sep=""))
	# Store list of fits and return params
	# params$list_fits <- list(abund=list_abund, survA=list_survA, survS=list_survL, catchA=list_catchA,
	#                          catchS=list_catchL, tagrecap=list_tagrecap)
	# return(params)
}
