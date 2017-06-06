#### Run CASAL Assessment

#' Run CASAL point estimate
#'
#' Run CASAL point estimate
#'
#' CASAL Commands
#'
#' -e: Calculate the point estimate
#'
#' -r: Run the population section once only and calculate the objective function
#'
#' -O: Output the estimated parameter values to the following files
#'
#' -f: Use a prefix on the names of the three input parameter files
#'
#' -p: Calculate likelihood or posterior profiles
#'
#' -q: Run quietly, suppress printing from within the population section.
#'
#' -i: Input one or more sets of free parameter values from a text file.
#' 	  With -r, run the model with each
#'
#' -P: Calculate projected outputs. Must use -i with name of file containing free
#' parameters. Results are dumped into outfile
#'
#' -g: Seed the random number generator with this positive integer value
#'
#' -s: [number prefix] Generate simulated observations
#'
#' @inheritParams casal
#' @param inputprefix prefix
#' @export
run_casal_pointest <- function(casal_path, mpd_dat, inputprefix, output_log, linux=0, intern = TRUE) {
	# CASAL can be called from anywhere if entry in:
	# System - Environmental variables:     Path: C:\Program Files\CASAL
	# Therefore, 'casal...' is enough as command here
	# casal -e -O mpd.dat -f prefix > output.log
	if(linux == 1) { 	# Linux
	  # The '2> /dev/null' is to re-direct stderr (casal output from the cmd line)
	  # so that it is not printed into the Condor out file
		run_casal <- paste("casal -e -Q -O ",casal_path,mpd_dat," -f ",casal_path,inputprefix,
		                   " > ",casal_path,output_log," 2> /dev/null", sep="")
    # print to screen
		print(run_casal)
		# run CASAL
		system(run_casal)
	} else {		# Windows
		run_casal <- paste("casal -e  -O ", casal_path,mpd_dat, " -f ",casal_path, inputprefix,
		                   " > ",casal_path,output_log, sep="")
		# print to screen
		print(run_casal)
		# run CASAL
		shell(run_casal, intern = intern)
	}
}

#' Run CASAL mvn
#'
#' Draw MVN samples for parameter estimates from CASAL model fit
#'
#' Based on code written by Steve Candy and used in WG-FSA-11/24
#' @inheritParams casal
#' @param nsamples number of samples (default=1000)
#' @export
run_casal_mvnsamples <- function(casal_path, mpd_dat, output_log, mvnsamples_dat, nsamples=1000)  {
  # mpd_mean is "mpd.dat"
  mpd_dat1 <- paste(casal_path,mpd_dat,sep="")
	# vcov is covariance matrix
	res_covmatrix <- casal::extract.covariance.matrix(file=output_log, path=casal_path)
	## MPD data
	mpd_mean <- matrix(scan(mpd_dat1,skip=1),ncol=length(scan(mpd_dat1,skip=1)),byrow=T)
	## Covariance matrix
	vcov <- as.matrix(res_covmatrix)
	# Make it robust to minor differences (lack of symmetry as required in rmvnorm)
	for (i in 1:nrow(vcov)) {
		for (j in 1:ncol(vcov)) {
			vcov[i,j] <- vcov[j,i]
		}
	}
	##!! ###***
	#isSymmetric(vcov)
	#vcov.se <- rep(0,length(mpd_mean))
	#for (r in 1:length(mpd_mean)){
	#	vcov.se [r] <- vcov [r,r] ^ 0.5
	#}
	## Create multivariate normal samples and store in CASALout_mvn_samples.dat
	mpd_mvn <- mvtnorm::rmvnorm(n=nsamples, mean=mpd_mean, sigma=vcov, method=c("svd"))
	mpd_mvn[mpd_mvn<0.001] <- 0.001      	# Protect against very low or negative values
	## Store mpd_mean and mvn results in mvn_samples.dat
	filename <- paste(casal_path, mvnsamples_dat, sep="")
	cat(scan(file=mpd_dat1, what=character(), nlines=1), file=filename)
	cat("\n",file=filename, append=TRUE)
	write.table(rbind(mpd_mean,mpd_mvn),file=filename, row.names=FALSE, col.names=FALSE, append=TRUE)
	##!! ###*** Selectivity confidence bounds
}

#' Run one set of projections using MVN samples
#'
#' Run one set of projections using MVN samples (using future_catches)
#' @inheritParams casal
#' @param proj_dat proj_dat???
#' @param set_seed set random number seed
#' @param rand_seed random number seed
#' @export
run_casal_proj_once <- function(casal_path, proj_dat, mvnsamples_dat, inputprefix="",
                                set_seed=0, rand_seed=0, linux=0) {
	# casal -P proj2000.dat -i mvn.samples.dat -f prefix
	if(set_seed == 1) {
		run_casal <- paste("casal -P ",casal_path,proj_dat," -Q -g ", rand_seed," -i ",
		                   casal_path,mvnsamples_dat," -f ",casal_path,inputprefix,sep="")
	} else {
		run_casal <- paste("casal -P ",casal_path,proj_dat," -Q -g 0 -i ",casal_path,
		                   mvnsamples_dat," -f ",casal_path,inputprefix,sep="")
	}
	if(linux == 1) { 	# Linux
	  # run CASAL. note "2> /dev/null" does not work here
		system(run_casal)
	} else {		# Windows
	  # run CASAL
	  shell(run_casal, intern = TRUE)
	}
}

#' Obtain likelihood profile for B0
#'
#' Obtain likelihood profile for B0
#' @inheritParams casal
#' @param nsteps number of MCMC steps
#' @param stepsize step size
#' @param est_file estimation file
#' @param profile_dat profile data
#' @param profile_out profile output file
#' @param profile_log profile log file
#' @export
run_casal_LLprofile <- function(nsteps=11, stepsize=10000, casal_path, inputprefix="", est_file,
                                output_log, profile_dat, profile_out, profile_log, linux=0) {
	## Change first @profile in estimation.csl to define appropriate lower and upper limits and number of steps (n, l, u)
	casalest <- casal::extract.csl.file(file=est_file, path=casal_path)
	res_quant <- casal::extract.quantities(file=output_log, path=casal_path)
	temp <- round(res_quant$'Scalar parameter values'$initialization.B0/stepsize,0)*stepsize
	casalest[["profile[1]"]]$n <- nsteps
	casalest[["profile[1]"]]$l <- temp - (nsteps-1)/2*stepsize
	casalest[["profile[1]"]]$u <- temp + (nsteps-1)/2*stepsize
	casal::write.csl.file(casalest,est_file,casal_path)
	## Create parameter estimates for each B0 (slow)
	# casal -p -q -O profile.dat > profile.out
	run_casal <- paste("casal -p -Q -O ",casal_path,profile_dat, " -f ",casal_path,
	                   inputprefix," > ",casal_path,profile_out, sep="")
	if(linux == 1) { 	# Linux
	  # run CASAL
		system(run_casal)
	} else {		# Windows
	  # run CASAL
	  shell(run_casal, intern = TRUE)
	}
	## Create LL estimates for each B0 based on parameters estimated in previous step (fast)
	# casal -i profile.dat -q -r > profile.log
	run_casal <- paste("casal -i ",casal_path,profile_dat, " -f ",casal_path,
	                   inputprefix," -Q -r > ",casal_path,profile_log, sep="")
	if(linux == 1) { 	# Linux
	  # run CASAL
		system(run_casal)
	} else {		# Windows
	  # run CASAL
	  shell(run_casal, intern = TRUE)
	}
}

#' Update CASAL data
#'
#' Update CASAL data
#' @inheritParams casal
#' @export
update_casal_data <- function(datass, Yr_current) {		# Setup function here, will be overwritten if running different scenarios
	# datass$sample_years <- array(0, dim=c(length(types),length(datass$years),length(datass$list_fishery)),
	#						dimnames=list(types=types,year=datass$years,fishery=datass$list_fishery))
	ycurr <- which(datass$years==Yr_current)
	age_y	<- which(datass$years==2000):ycurr
	surv_y <- seq(6,ycurr,2)
	datass$sample_years["catchlen", (ycurr-5):ycurr,] <- 0
	datass$sample_years["catchage", age_y,] <- 1
	datass$sample_years["survey", surv_y,2] <- 1
	datass$sample_years["surveyage", surv_y,] <- 1
	datass$sample_years["surveylen", (ycurr-10):ycurr,] <- 0
	datass$sample_years["cpue", (ycurr-10):ycurr,] <- 0	# Relative index (cpue)
	datass$sample_years["tagging", (ycurr-10):(ycurr-1),] <- 1	# No tags in last year
	return(datass)
}

#' Update future catches in population.csl
#'
#' Update future catches in population.csl
#' @inheritParams casal
#' @export
update_future_catches <- function(datass, TAC, Yr_current, use_specific_ratios=0,
                                  specific_ratios=1, casal_path="", pop_csl="population.csl") {
	## Get population file
	casalpop	<- casal::extract.csl.file(paste(casal_path, pop_csl, sep=""))
	## Calculate ratios
	if(use_specific_ratios == 1) {		# Use specific_ratios
		ratios <- as.numeric(specific_ratios)
	} else {							# Estimate ratios from catches in assessment year
		refy <- which(datass$years==Yr_current)
		catch_refy <- vector(mode="numeric", length=length(datass$list_fishery))
		count <- 0
		for (met in 1:length(datass$list_metier)) {
			for (ss in 1:datass$season[3]){
				for (rr in 1:datass$region[3]){
					if(apply(datass$catch[[met]][,refy,,ss,rr,],c(4:5),sum) > 0) {
						count <- count + 1
						catch_refy[count] <- as.vector(apply(datass$catch[[met]][,refy,,ss,rr,],c(2,4:5),sum))
					}
				}
			}
		}
		# Take ratios of reference year for projections
		ratios <- round(catch_refy / sum(catch_refy), 2)
	}
	if(sum(ratios) > 1) ratios <- ratios/sum(ratios)
	## Update future catches in population.csl
	proj_catch <- TAC * ratios
	for (i in 1:length(datass$list_fishery)) {
		fisheryN <- paste("fishery[",datass$list_fishery[i],"]",sep="")
		casalpop[[fisheryN]]$future_constant_catches <- proj_catch[i]
	}
	# write CASAL population file
	casal::write.csl.file(casalpop,pop_csl,casal_path)
	print("Future catches updated")
}

#' Run CASAL assessment (with adjustment for process error, with/witout TAC finder)
#'
#' Run CASAL assessment (with adjustment for process error, with/witout TAC finder)
#'
#' Function to run a complete CASAL assessment loop, including:
#'
#' 	- Create input data (source(file="M3_CASAL1_Data.r"))
#'
#'	- Create CASAL input files (create_casal_file_pop, create_casal_file_est,
#'	create_casal_file_out)
#'
#'	- Run CASAL model for point estimatation (run_casal_pointest)
#'
#'	- Draw MVN samples for parameter estimates (run_casal_mvnsamples)
#'
#'	- Find 'optimal' TAC for next year (run_casal_TAC_finder to update catches,
#'	run projections to select optimal TAC that satifies the CCAMLR decision rules)
#'
#' ctrl Parameters for general control
#'
#' om	Parameters for Operating model
#'
#' datass	Parameters for Assessmment
#'
#' Yr_current	Assessment year, also used to calculate TAC ratios for
#' all fisheries (if use_specific_ratios=0)
#'
#' TAC Initial catch level
#'
#' @inheritParams casal
#' @param ctrl Parameters for general control
#' @param om Parameters for Operating model
#' @param tag tag data/parameters
#' @param datass Parameters for Assessmment
#' @param TAC_init Initial catch level
#' @export
run_complete_casal_assessment <- function(ctrl, om, tag, datass, Yr_current, TAC_init=0, intern=TRUE) {
  # extract parameters
	use_specific_ratios <- ctrl$use_specific_ratios				# If 1, then use TAC ratios provided by 'specific_ratios'
	specific_ratios <- datass$Fish[,"ProjCatch"]			# Specific TAC ratios (if use_specific_ratios=1)
	casal_path <- ctrl$casal_path					    # File path
	inputprefix <- ctrl$inputprefix					    # Prefix of input files (population.csl, estimation.csl and output.csl files)
	outputprefix <- ctrl$outputprefix				    # Prefix of output files
	pop_csl <- ctrl$pop_csl							# File name for population.csl
	est_csl <- ctrl$est_csl							# File name for estimation.csl
	out_csl <- ctrl$out_csl							# File name for output.csl
	output_log <- ctrl$output_log					    # File name for output.log
	pop_skel_csl <- ctrl$pop_skel_csl					# File name for population_skel.csl
	est_skel_csl <- ctrl$est_skel_csl					# File name for estimation_skel.csl
	mpd_dat <- ctrl$mpd_dat						    # File name for mpd.dat
	proj_dat <- ctrl$proj_dat					    # File name for projection.dat
	mvnsamples_dat <- ctrl$mvnsamples_dat				    # File name for input mvn samples
	linux <- ctrl$linux							# Run on linux (1) or Windows (0)
	pin_TAC_finder <- ctrl$pin_TAC_finder			        # Find TAC (1=Yes, 0=No)
	ssb_level_deplet <- ctrl$ssb_level_deplet			    # Reference SSB level for depletion (0.2)
	ssb_level_target <- ctrl$ssb_level_target			    # Reference SSB level for target (e.g. 0.5 for toothfish, 0.75 for icefish)
	ref_prob_deplet <- ctrl$ref_prob_deplet				    # Refernce probability of SSB dropping below ssb_level_deplet (0.1)
	ref_prob_target <- ctrl$ref_prob_target				    # Reference measure of SSB that needs to be at/above ssb_level_target ('median', not used)
	tol <- ctrl$TAC_Finder_tol					# Tolerance level defining the acceptable target reference levels for the search,
																#    Target range is from (ssb_level_target) to (ssb_level_target + tol)
																#    Depletion range is from (ref_prob_deplet - tol) to (ref_prob_deplet)
	TAC_multiplier <- ctrl$TAC_Finder_TAC_multiplier		# Provides multiplier for in/decrementing TAC values (e.g. 1.2)
	iterations_main <- ctrl$TAC_Finder_iterations_main		# Number of iterations when 'zooming in' on specific TAC level
	set_seed <- ctrl$TAC_Finder_set_seed				# If = 1, seed the random number generator with rand_seed
	rand_seed <- ctrl$TAC_Finder_rand_seed			# Positive integer value to seed random number generator if set_seed=1
	#print(dir())
	## !! ###*** consider adding file.exists()
	if(file.exists(paste0(casal_path,mpd_dat))) file.remove(paste0(casal_path,mpd_dat))
	if(file.exists(paste0(casal_path,output_log))) file.remove(paste0(casal_path,output_log))
	## Run CASAL model for point estimatation
	print("Run Casal for point estimate")
	assesstime <- proc.time()		# For time measurement
	###*** can we remove this?
	# use tryCatch to avoid 'Execution halted' when CASAL cannot find a solution
	# return(tryCatch(run_casal_pointest(casal_path, inputprefix=inputprefix, mpd_dat=mpd_dat, output_log=output_log,
	#	       error=function(e) print("CASAL failed to find solution")))
  # tryCatch(run_casal_pointest(casal_path, inputprefix=inputprefix, mpd_dat=mpd_dat, output_log=output_log, linux=linux),
  #   		   error=function(e) print("CASAL failed to find solution"))
	tryCatch(run_casal_pointest(casal_path, inputprefix=inputprefix, mpd_dat=mpd_dat,
	                            output_log=output_log, linux=linux, intern=intern),
	         error=function(e) print("CASAL failed to find solution"))
	#run_casal_pointest(casal_path, inputprefix=inputprefix, mpd_dat=mpd_dat, output_log=output_log)
	print("Total time for CASAL Assessment: ")
	print(proc.time() - assesstime)		# For time measurement
  ## Adjust for process error (find cv for process error, adjust Effective Sample Size, rerun CASAL)
	print("Adjust for process error not needed for these scenarios")
	# If results found, then mpd_dat is created if it does not exist and filled with some data
	# It is faster to check mpd_dat than output_log
	## Save output_log by year of assessment
	nname <- unlist(strsplit(output_log,"\\."))
	file.copy(from=paste(casal_path,output_log,sep=""),
	          to=paste(casal_path,nname[1],Yr_current,".",nname[2],sep=""),
	          overwrite = TRUE)
  ## calculate TAC
	if(pin_TAC_finder== 1) {
		if(file.exists(paste(casal_path,mpd_dat,sep=""))) {
			if(length(scan(file=paste(casal_path,mpd_dat,sep=""),what="character",skip=1)) > 0){
			  ## Draw MVN samples for parameter estimates
				print("Run mvn sampling")
				run_casal_mvnsamples(casal_path=casal_path, mpd_dat=mpd_dat, output_log=output_log,
									 mvnsamples_dat=mvnsamples_dat, nsamples=datass$output[["n_projections"]])
				## Find 'optimal' TAC for next year
				print(paste("Run TAC finder for ",Yr_current+1,sep=""))
				res <- run_casal_TAC_finder(datass=datass, TAC=TAC_init, Yr_current, use_specific_ratios, specific_ratios,
					casal_path=casal_path, inputprefix, pop_csl, proj_dat, mvnsamples_dat, linux,
					ssb_level_deplet, ssb_level_target, ref_prob_deplet, ref_prob_target,
					tol, TAC_multiplier, iterations_main, set_seed, rand_seed)
				# Return TAC for next year
				return(res$Final_TAC)
			} else {
			  # If CASAL fit was unsuccessful, return same TAC as in previous year
				return(TAC_init)
			}
		} else {
		  # If CASAL fit was unsuccessful, return same TAC as in previous year
			return(TAC_init)
		}
	} else {
	  # If pin_TAC_finder==0, return NA: TAC will then not be used in subsequent year as catch level
		return(NA)
	}
}

#' CASAL TAC Finder
#'
#' CASAL TAC Finder
#' @inheritParams casal
#' @param proj_dat project data
#' @param TAC_multiplier TAC multiplier
#' @param iterations_main iterations
#' @export
run_casal_TAC_finder <- function(datass, TAC, Yr_current, use_specific_ratios=0, specific_ratios=1,
                                 casal_path, inputprefix, pop_csl, proj_dat, mvnsamples_dat, linux,
                                 ssb_level_deplet=0.2, ssb_level_target=0.5, ref_prob_deplet=0.1,
                                 ref_prob_target="median", tol = 0.01, TAC_multiplier=1.2, iterations_main=20,
                                 set_seed=0, rand_seed=0) {
	##!! ###*** To do:  Simplify
	# Better alternative routine than 'while?
	# Search function needs to be tested for robustness
	print("Start with initial TAC")
	## Search with extrapolation only
	## 1. Use starting TAC and extrapolate from initial TAC value to get quickly to the vicinity of the 'optimal' TAC
	counter <- 0
	end_extrapolation <- FALSE
	while(end_extrapolation == FALSE) {
		counter <- counter + 1
		res <- run_casal_projections(datass, TAC, Yr_current, use_specific_ratios,
		                             specific_ratios, casal_path, inputprefix, pop_csl,
		                             proj_dat, mvnsamples_dat, linux, ssb_level_deplet,
		                             ssb_level_target, ref_prob_deplet, ref_prob_target,
		                             tol, set_seed, rand_seed)
		print(paste("ProbTest: ",res$Prob_test, "; TAC: ",round(TAC,0),";   P(Depletion) ",
		            round(res$Probabilities["Depletion"],3), ";   P(Target) ",
		            round(res$Probabilities["Target"],3), ";   ExtrapolCatch: ",
		            round(min(res$Extrapolated_Catch),0),sep=""))
		## If solution found:
		if(res$Prob_test == 1) break			# Solution found (at least one prob within tol range, other probability ok)
		## If some problem:
		if(res$Prob_test == 5) {				# Stop search (some error) and use slightly smaller TAC
			TAC <- 0.9 * TAC
			print("Some error occurred - use slightly smaller TAC in next step")
			break
		}
		## If TAC should be smaller, but evaluated TAC is already very small:
		if(res$Prob_test != 2 & TAC <= 10) {		# Stop search (both probs are bad and TAC is already close to 0)
			TAC <- 10		 # Has to be >0 otherwise no data collection
			print("Only a minimum catch of 10t can be taken")
			break
		}
		## If TAC should be either smaller or larger: Use min extrapolated catch in next step
		oldTAC <- TAC
		TAC <- max(round(min(res$Extrapolated_Catch),0),0) 	# Take lower of extrapolated estimates for new TAC, but at least 0
		# If TAC = NaN (with oldTAC = 0, extrapolation returns NaN), take arbitrary value of 1000
		if (TAC == "NaN") TAC <- 1000
		# If TAC is very large or Inf, use arbitrary new TAC of 5-10*oldTAC
		if (TAC >= 10*oldTAC | TAC == Inf) TAC <- runif(1,min=5,max=10)*oldTAC
		# If both probs are 'good' but TAC is 0 (e.g. when extrapolor catch was negative with depletion prob = 0 and target prob > 1)
		if (res$Prob_test == 2 & TAC == 0) TAC <- runif(1,min=5,max=10)*oldTAC		# Increase TAC
		## Rules for termination
		if(counter < iterations_main) {
			print("Use TAC from extrapolation in next step")
		}
		if (counter == iterations_main) {	# Sometimes search get stuck in a loop - change TAC randomly and allow for another 10 loops
			print("Number of main iterations reached: Change TAC randomly & allow 10 more iterations")
			TAC <- round(TAC * abs(rnorm(1, mean=1, sd=0.6)),0)
		}
		if(counter > iterations_main + 10) {
			print("Too many iterations - process ended")
			break
		}
	}
	list(Final_TAC=round(TAC,0), Probabilities=res$Probabilities)
}

#' Run CASAL projections
#'
#' Update catches, run projections once, return probabilities & extrapolated catch
#' and check whether probabilities are within tolerance range
#' @inheritParams casal
#' @param use_specific_ratios switch for specfic ratios
#' @param specific_ratios specific ratios
#' @param proj_dat projection data?
#' @export
run_casal_projections <- function(datass, TAC, Yr_current, use_specific_ratios=0, specific_ratios=1,
                                  casal_path, inputprefix, pop_csl, proj_dat, mvnsamples_dat, linux,
                                  ssb_level_deplet=0.2, ssb_level_target=0.5, ref_prob_deplet=0.1,
                                  ref_prob_target="median", tol = 0.01, set_seed=0, rand_seed=0) {
  ## Update population.csl with new future_catches
	update_future_catches(datass, TAC, Yr_current, use_specific_ratios,
	                      specific_ratios, casal_path, pop_csl)
	## Run projections
	run_casal_proj_once(casal_path, proj_dat, mvnsamples_dat,
	                    inputprefix, set_seed, rand_seed, linux)
	## Calculate depletion and target probabilities
	print("Calculate depletion and target probabilities")
	col_names <- scan(paste(casal_path, proj_dat, sep=""),what="character",skip=8,nlines=1)
	#ssb_data <- read.table(paste(casal_path,proj_dat,sep=""),skip=9)
	tryCatch(ssb_data	<- read.table(paste(casal_path, proj_dat, sep=""), skip=9),
								error=function(e) print("Error in proj.dat file - no calculation of depletion or target probabilities"))
	if (exists("ssb_data")) {
		ssb_data <- ssb_data[,which(substring(col_names,1,3)=="SSB")]
		colnames(ssb_data) <- col_names[which(substring(col_names,1,3)=="SSB")]
		res <- get_ssb_probabilities(abs_data=ssb_data, datass_yr=datass$year, Catch=TAC,
		                             ssb_level_deplet, ssb_level_target, ref_prob_deplet,
		                             ref_prob_target)
		## Testing for depletion and target probabilities
		obs_prob_deplet <- round(res$Probabilities["Depletion"],3)
		obs_prob_target <- round(res$Probabilities["Target"],3)
		# Endpoints found if obs_prob_deplet is within tolerance levels and obs_prob_deplet is above reference level, or
		#                 if obs_prob_target is within tolerance levels and obs_prob_deplet is below reference level
		tol_deplet <- c(ref_prob_deplet - tol, ref_prob_deplet)		# Reference: Probability (e.g. 0.1)
		tol_target <- c(ssb_level_target - 0.5*tol, ssb_level_target + tol)	# Reference: Median being at ssb_level_target (e.g. 0.5)
		## Tests:
		if (is.na(obs_prob_deplet) | is.na(obs_prob_target)) {
			print("Error in calculation of depletion or target probabilities")
			list(Prob_test=5, Probabilities=0, Extrapolated_Catch=0)	# Error
		} else {
			if (obs_prob_deplet < tol_deplet[1] & obs_prob_target > tol_target[2]) prob_test <- 2		# Both prob are 'too good': increase catch
			if (obs_prob_deplet > tol_deplet[2] & obs_prob_target < tol_target[1]) prob_test <- 4		# Both prob are 'bad': decrease catch
			if (obs_prob_deplet > tol_deplet[2] | obs_prob_target < tol_target[1]) prob_test <- 3		# At least one prob is 'bad'
			if (obs_prob_deplet >= tol_deplet[1] & obs_prob_deplet <= tol_deplet[2] & obs_prob_target >= tol_target[1]) prob_test <- 1
			if (obs_prob_target >= tol_target[1] & obs_prob_target <= tol_target[2] & obs_prob_deplet <= tol_deplet[2]) prob_test <- 1
			# Return values
			list(Prob_test=prob_test, Probabilities=res$Probabilities, Extrapolated_Catch=res$Extrapolated_Catch)
		}
	} else {	# If error further above, return prob_test = 5
		list(Prob_test=5, Probabilities=0, Extrapolated_Catch=0)
	}
}

#' Get projections for absolute and relative SSB
#'
#' Get projections for absolute and relative SSB.
#' @param abs_data Input matrix with rows = trials & cols = years,
#' data = absolute data to be plotted as box plots
#' @param datass_yr	datass$year
#' @param Catch Catch level (used in extrapolation)
#' @param ssb_level_deplet Reference SSB level for depletion (e.g. 0.2)
#' @param ssb_level_target Reference SSB level for target (e.g. 0.5 for
#' toothfish)
#' @param ref_prob_deplet	Reference probability of SSB dropping below
#' ssb_level_deplet
#' @param ref_prob_target	Reference measure of SSB that needs to be at/above
#' ssb_level_target (median). Using the median is hardwired, ref_prob_target
#' is not actually used.
#' @export
get_ssb_probabilities <- function(abs_data, datass_yr, Catch, ssb_level_deplet,
                                  ssb_level_target, ref_prob_deplet, ref_prob_target) {
	# Relative SSB
	init <- abs_data[,1]		# Note: these are not the B0s in the MVN file
	rel_data <- sweep(abs_data,c(1), init,"/")
	# Return probabilities & extrapolated catch that would satisfy decision rule
	# assuming that zero catch has depletion prob=0 and target prob=1)
	rel_data_proj <- rel_data[,(datass_yr[3]+1):ncol(rel_data)]	# Projection years only
	obs_prob_deplet <- sum(rowSums(rel_data_proj <= ssb_level_deplet)>0)/nrow(rel_data)
	obs_prob_target <- median(rel_data_proj[,ncol(rel_data_proj)])		# 'median' is ref_prob_target
	extrapol_catch <- round(c((ref_prob_deplet  - 0)/(obs_prob_deplet-0)*(Catch-0),
					           (ssb_level_target - 1)/(obs_prob_target-1)*(Catch-0)),0)
	list(Probabilities = c(Depletion=obs_prob_deplet, Target=obs_prob_target),
		 Extrapolated_Catch = c(Depletion_rule=extrapol_catch[1], Target_rule=extrapol_catch[2]))
}
