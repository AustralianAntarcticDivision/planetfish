##### Create population.csl file for CASAL Assessment

#' Create CASAL population.csl file
#'
#' Create CASAL population.csl file
#' @inheritParams casal
#' @export
create_casal_file_pop <- function(params, casal_path, skel_csl, csl) {
  #### This funtion sets up a casal.population file for an assessment
  # params: Object with input parameters for CASAL assessment para[["ass"]]
	# 		  Data includes catch and fish tagged
	# params    	<- datass
	# casal_path 	<- ""
	# skel_csl  	<- "casal_population_skel.csl"
	# csl			<- "casal_population.csl"
	#### Steps of this function:
	# extract.csl.file(xxx)							# 0. Read in csl skeleton file
	# casalpop[["test1"]] <- NULL					# 2. Remove list elements
	# casalpop <- c(casalpop,"test3"=list(a))		# 3. Create new list elements (steps 1-3 is to maintin order of list elements
													#    even if when more than one fishery)
	# casalpop <- c(casalpop,"test2"=list(casalpop[["test1"]]))	# or: duplicate list elements
													# 4. Fill in data to list elements
  # extract casal parameters from skeleton file
	casalpop <- casal::extract.csl.file(paste(casal_path, skel_csl, sep=""))
	## Get assessment year
	ycurr <- which(params$years==params$Yr_current)
	#### 0. Prepare casalpop file (for sequence, see Manual Chapter 14.1.2 The input parameter files)
	# Remove some list elements (that relate to individual fisheries or tagging data)
	casalpop[["stock_names"]] <- NULL
	casalpop[["fishery[fishery1]"]] <- NULL
	casalpop[["selectivity_names"]] <- NULL
	casalpop[["selectivity[fishery1]"]] <- NULL
	casalpop[["n_tags"]] <- NULL
	casalpop[["tag_names"]] <- NULL
	casalpop[["tag_shedding_rate"]] <- NULL
	casalpop[["tag_loss_props"]] <- NULL
	casalpop[["tag_growth_loss[2000Tags]"]] <- NULL
	casalpop[["tag[2000Tags]"]] <- NULL
	#### 1. Initialization
	casalpop[["initialization"]]$B0 <- as.integer(params$B0)
	#### 2. Partition
	casalpop[["size_based"]]$value <- params$size_based
	casalpop[["min_age"]]$value <- params$age[1]
	casalpop[["max_age"]]$value <- params$age[2]
	casalpop[["plus_group"]]$value <- params$age_plus_group
	casalpop[["sex_partition"]]$value <- params$sex_partition
	casalpop[["mature_partition"]]$value <- params$mature_partition
	casalpop[["n_areas"]]$value <- length(params$regions)
	casalpop[["area_names"]]$value <- c(params$regions)
	casalpop[["n_stocks"]]$value <- params$n_stocks
	if (params$n_stocks > 1) {
		casalpop[["stock_names"]]$command <- "stock_names"
		casalpop[["stock_names"]]$value <- params$stock_names
	}
	#### 3. Time sequence
	casalpop[["initial"]]$value <- params$year[1]
	casalpop[["current"]]$value <- params$year[2]
	casalpop[["final"]]$value <- params$year[2] + params$y_proj
	casalpop[["annual_cycle"]]$command <- "annual_cycle"
	casalpop[["annual_cycle"]]$value <- character(0)
	casalpop[["annual_cycle"]]$time_steps <- params$time_steps
	casalpop[["annual_cycle"]]$recruitment_time <- params$recruitment_time
	casalpop[["annual_cycle"]]$recruitment_areas <- params$recruitment_areas[1]	# Only 1 area per stock
	casalpop[["annual_cycle"]]$spawning_time <- params$spawning_time
	casalpop[["annual_cycle"]]$spawning_part_mort <- params$spawning_part_mort
	if (params$SSB_for_all_areas_pin == 1) {
		casalpop[["annual_cycle"]]$spawning_all_areas <- params$spawning_all_areas
	} else {
		casalpop[["annual_cycle"]]$spawning_areas <- params$spawning_areas
	}
	casalpop[["annual_cycle"]]$spawning_p <- params$spawning_p
	casalpop[["annual_cycle"]]$spawning_use_total_B <- params$spawning_use_total_B
	casalpop[["annual_cycle"]]$aging_time <- params$aging_time
	casalpop[["annual_cycle"]]$growth_props <- params$growth_props
	casalpop[["annual_cycle"]]$M_props <- params$M_props
	casalpop[["annual_cycle"]]$baranov <- params$baranov
	casalpop[["annual_cycle"]]$fishery_names <- params$list_fishery
	casalpop[["annual_cycle"]]$fishery_times <- params$list_season
	casalpop[["annual_cycle"]]$fishery_areas <- params$list_region
	casalpop[["annual_cycle"]]$n_migrations <- params$n_migrations
	if (length(params$n_migrations) > 0) {
		casalpop[["annual_cycle"]]$migration_names <- params$migration_names
		casalpop[["annual_cycle"]]$migration_times <- params$migration_times
		casalpop[["annual_cycle"]]$migrate_from <- params$migrate_from
		casalpop[["annual_cycle"]]$migrate_to <- params$migrate_to
	}
	#### 4. Recruitment
	casalpop[["y_enter"]]$value <- params$rec_y_enter
	casalpop[["standardise_YCS"]]$value <- params$rec_standardise_YCS
	casalpop[["recruitment"]]$command <- "recruitment"
	casalpop[["recruitment"]]$value <- character(0)
	casalpop[["recruitment"]]$YCS_years <- params$rec_YCS_years
	casalpop[["recruitment"]]$YCS <- params$rec_YCS
	casalpop[["recruitment"]]$first_free <- params$rec_first_free
	casalpop[["recruitment"]]$last_free <- params$rec_last_free
	casalpop[["recruitment"]]$year_range <- params$rec_year_range
	casalpop[["recruitment"]]$SR <- params$rec_SR
	casalpop[["recruitment"]]$sigma_r <- params$rec_sigma
	casalpop[["recruitment"]]$steepness <- params$rec_steepness
	casalpop[["recruitment"]]$rho <- params$rec_rho
	casalpop[["recruitment"]]$p_male <- params$rec_p_male
	#### 5. Recruitment variability
	casalpop[["randomisation_method"]]$value <- params$rec_randomisation_method
	casalpop[["first_random_year"]]$value <- params$rec_first_random_year
	#### 6. Maturation
	casalpop[["maturity_props"]]$command <- "maturity_props"
	casalpop[["maturity_props"]]$value <- character(0)
	if (params$by_sex==0) {
		casalpop[["maturity_props"]]$all <- params$maturity_props_all
	} else {
		asses[["maturity_props"]]$female <- params$maturity_props_female
		asses[["maturity_props"]]$male <- params$maturity_props_male
	}
	#### 7. Migration
	if (params$n_migrations > 0) {
		for (ii in 1:params$n_migrations) {
			migrN <- paste("migration[",params$migration_names[ii],"]",sep="")
			casalpop[[migrN]]$command <- paste("migration ",params$migration_names[ii],sep="")
			casalpop[[migrN]]$value <- character(0)
			casalpop[[migrN]]$migrators <- params$migrators[ii]
			casalpop[[migrN]]$rates_all	<- params$rates_all[[ii]]
		}
	}
	#### 8. Population file: Natural mortality
	casalpop[["natural_mortality"]]$command <- "natural_mortality"
	casalpop[["natural_mortality"]]$value <- character(0)
	if (params$by_sex==0) {
		casalpop[["natural_mortality"]]$all <- params$estnatM[[1]][1]
	} else {
		casalpop[["natural_mortality"]]$female <- params$estnatM[[1]][1]
		casalpop[["natural_mortality"]]$male <- params$estnatM[[2]][2]
	}
	#### 9. Fishing / catch
	for (ff in 1:length(params$list_fishery)) {
		fishery <- params$list_fishery[ff]
		if (sum(params$catch[[fishery]]) > 0) {
			fisheryN <- paste("fishery[",fishery,"]",sep="")
			casalpop[[fisheryN]]$command <- paste("fishery ",fishery,sep="")
			casalpop[[fisheryN]]$value <- character(0)
			casalpop[[fisheryN]]$years <- dimnames(params$catch[[fishery]])$year
			casalpop[[fisheryN]]$catches <- as.vector(apply(params$catch[[fishery]],c(2),sum))
			casalpop[[fisheryN]]$U_max <- params$U_max
			casalpop[[fisheryN]]$selectivity <- params$list_sel[ff]
			casalpop[[fisheryN]]$future_constant_catches <- params$future_constant_catches[ff]
			#casalpop[[fisheryN]]$future_years			<-
			#casalpop[[fisheryN]]$future_catches		<-
		}
	}
	#### 10. Selectivities
	casalpop[["selectivity_names"]]$command <- "selectivity_names"
	casalpop[["selectivity_names"]]$value <- params$selectivity_names
	for (ii in 1:length(params$selectivity_names)) {
		selN <- paste("selectivity[",params$selectivity_names[ii],"]",sep="")
		casalpop[[selN]]$command <- paste("selectivity ",params$selectivity_names[ii],sep="")
		casalpop[[selN]]$value <- character(0)
		if (params$by_sex==0) {
			casalpop[[selN]]$all <- params$selN_all[[ii]]
		} else {
			casalpop[[selN]]$female <- params$selN_female[[ii]]
			casalpop[[selN]]$male <- params$selN_male[[ii]]
		}
	}
	#### 11. Size at age
	casalpop[["size_at_age_type"]]$value <- params$size_at_age_type
	casalpop[["size_at_age_dist"]]$value <- params$size_at_age_dist
	casalpop[["size_at_age"]]$command <- "size_at_age"
	casalpop[["size_at_age"]]$value <- character(0)
	if (params$by_sex==0) {
		casalpop[["size_at_age"]]$Linf <- params$estgrowth[[1]][1]		# if sex partition: Linf_male or Linf_female etc
		casalpop[["size_at_age"]]$k <- params$estgrowth[[1]][2]
		casalpop[["size_at_age"]]$t0 <- params$estgrowth[[1]][3]
		casalpop[["size_at_age"]]$cv <- params$estgrowth[[1]][4]
	} else {
		casalpop[["size_at_age"]]$Linf_female <- params$estgrowth[[1]][1]
		casalpop[["size_at_age"]]$k_female <- params$estgrowth[[1]][2]
		casalpop[["size_at_age"]]$t0_female <- params$estgrowth[[1]][3]
		casalpop[["size_at_age"]]$cv_female <- params$estgrowth[[1]][4]
		casalpop[["size_at_age"]]$Linf_male <- params$estgrowth[[2]][1]		# if sex partition: Linf_male or Linf_female etc
		casalpop[["size_at_age"]]$k_male <- params$estgrowth[[2]][2]
		casalpop[["size_at_age"]]$t0_male <- params$estgrowth[[2]][3]
		casalpop[["size_at_age"]]$cv_male <- params$estgrowth[[2]][4]
	}
	#### 12. Size-weight
	casalpop[["size_weight"]]$command <- "size_weight"
	casalpop[["size_weight"]]$value <- character(0)
	if (params$by_sex==0) {
		casalpop[["size_weight"]]$a <- params$estWL[[1]][1] 		# if sex partition: a_male or a_female etc
		casalpop[["size_weight"]]$b <- params$estWL[[1]][2]
	} else {
		casalpop[["size_weight"]]$a_female <- params$estWL[[1]][1]
		casalpop[["size_weight"]]$b_female <- params$estWL[[1]][2]
		casalpop[["size_weight"]]$a_male <- params$estWL[[2]][1]
		casalpop[["size_weight"]]$b_male <- params$estWL[[2]][2]
	}
	casalpop[["size_weight"]]$verify_size_weight <- params$verify_size_weight
	#### 13. Population file: Tag release events
	## Tag release data in params$tag_props_all and params$tag_numbers is already summarised by assessment fisheries & regions
	## params$tag_names: Stores tag partitions & release events
	## TagName: Tag partition: Tag releases in particular year & region, e.g. Tags2009_R1, link releases with recaptures
	## TagEvent: Tag partition released by a particular fishery (and its selectivity)
	ttag <- params$sample_years["tagging",1:ycurr,,drop=FALSE]
	if (sum(ttag) > 0) {
		tag_names <- params$tag_names
		## General tag information
		casalpop[["n_tags"]]$command <- "n_tags"
		casalpop[["n_tags"]]$value <- nrow(tag_names)
		casalpop[["tag_names"]]$command <- "tag_names"
		casalpop[["tag_names"]]$value <- tag_names[,"TagName"]
		casalpop[["tag_shedding_rate"]]$command <- "tag_shedding_rate"
		casalpop[["tag_shedding_rate"]]$value <- rep(params$tag_shedding_rate,nrow(tag_names))
		casalpop[["tag_loss_props"]]$command <- "tag_loss_props"
		casalpop[["tag_loss_props"]]$value <- params$tag_loss_props
		for (ii in 1:nrow(tag_names)) {
			tagN <- paste("tag_growth_loss[",tag_names[ii,"TagName"],"]",sep="")
			casalpop[[tagN]]$command <- paste("tag_growth_loss ",tag_names[ii,"TagName"],sep="")
			casalpop[[tagN]]$value <- character(0)
			casalpop[[tagN]]$nogrowth_period <- params$nogrowth_period
		}
		## Tag-release events
		for (ff in 1:nrow(tag_names)) {
			tevent <- tag_names[ff,"TagEvent"]
			tyear	<- tag_names[ff,"Year"]
			tfish	<- tag_names[ff,"Fishery"]
			tagN <- paste("tag[",tevent,"]",sep="")
			# Only tag releases if there has been some catch taken that year
			if(sum(params$catch[[tfish]][,tyear,,]) > 0) {
				props <- params$tag_props_all[[tfish]][,tyear,,,,drop=FALSE]
				if (params$by_sex == 0)  props	<- apply(props,c(1,2),sum)
				# Store in casalpop
				casalpop[[tagN]]$command <- paste("tag ",tevent,sep="")
				casalpop[[tagN]]$value <- character(0)
				casalpop[[tagN]]$tag_name <- tag_names[ff,"TagName"]
				casalpop[[tagN]]$release_type <- params$tag_release_type
				casalpop[[tagN]]$area <- as.vector(tag_names[ff,"Region"])
				#casalpop[[tagN]]$stock 			<-
				casalpop[[tagN]]$sex <- params$tag_sex
				casalpop[[tagN]]$year <- tyear
				casalpop[[tagN]]$step <- params$Fish[params$Fish[,"Fishery"]==tfish,"Season"]
				casalpop[[tagN]]$mature_only <- params$tag_mature_only
				casalpop[[tagN]]$ogive <- params$Fish[params$Fish[,"Fishery"]==tfish,"Sel"]
				casalpop[[tagN]]$number <- as.vector(round(params$tag_numbers[[tfish]][tyear],0))
				casalpop[[tagN]]$plus_group <- params$len_plus_group
				casalpop[[tagN]]$class_mins <- params$class_mins
				if (params$by_sex == 0) {
					casalpop[[tagN]]$props_all <- paste(round(props,4),collapse=" ")
				} else {
					casalpop[[tagN]]$props_female <- paste(round(props[,,dimnames(props)$sex=="f",,],4),collapse=" ")
					casalpop[[tagN]]$props_male <- paste(round(props[,,dimnames(props)$sex=="m",,],4),collapse=" ")
				}
				casalpop[[tagN]]$mortality <- params$tag_mortality
			}
		}
	}
	# write the casal population file
	casal::write.csl.file(casalpop,paste(casal_path,csl,sep=""))
}
