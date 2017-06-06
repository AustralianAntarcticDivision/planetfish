## functions for sampling the fished population, age, length structures, surveys
## etc (moved from M2_Functions.R)

#' Sampling
#'
#' Sampling function
#' @param type sampling type
#' @param om operating model object
#' @param pop population object
#' @param sampling sampling object
#' @param fleet fleet object
#' @param obs observations?
#' @param ff fleets???
#' @param y year???
#' @param ss ss???
#' @param rr region???
#' @export
sample_survey	<- function(type, om, pop, sampling, fleet, obs, ff, y, ss, rr) {
  ## Sampling of catch-at-age
  ages1 <- om$ages + om$growth_props[ss]
  lenbins <- sampling$len_classes
  ageing_cv <- sampling$ageing_error
  # Storage
  surveyage <- obs[[ff]]$survey_age_n[,y,,ss,rr] * 0
  surveylen <- obs[[ff]]$survey_len_n[,y,,ss,rr] * 0
  survey <- obs[[ff]]$survey_index[,y,,ss,rr]
  for (i in 1:length(om$sex)) {
    ## Survey Index: Population.n * fishery selectivity * Survey catchability * weight, adjusted for overall variability
    agecomp <- pop$n[,y,i,ss,rr] * fleet[[ff]]$landings_sel[,y,i,ss,rr] * sampling$survey_q  # Numbers at age available
    survey[i] <- sum(agecomp * pop$wt[,y,i,ss,rr])
    if (sampling$survey_var_units == "normal") survey[i] <- survey[i] * (rnorm(1,mean=1,sd=sampling$survey_var))
    if (sampling$survey_var_units == "lognormal") survey[i] <- survey[i] * exp(rnorm(1,mean=0,sd=sampling$survey_var))
    ## Surveys at length/age: Adjust agecomp for variability with cvs for
    #	(1) overall numbers (the same for all ages/lengths)
    # (2) individual numbers (by age/length)
    # (1) Adjust for overall numbers
    # if (sampling$survey_var_units == "normal")    agecomp <- agecomp *    (rnorm(1,mean=1,sd=sampling$survey_var_N))
    # if (sampling$survey_var_units == "lognormal") agecomp <- agecomp * exp(rnorm(1,mean=0,sd=sampling$survey_var_N))
    # (2) Adjust for individual numbers
    if (sampling$survey_var_units == "normal") agecomp <- agecomp * (rnorm(length(agecomp),mean=1,sd=sampling$survey_var_n))
    if (sampling$survey_var_units == "lognormal") agecomp <- agecomp * exp(rnorm(length(agecomp),mean=0,sd=sampling$survey_var_n))
    if(ageing_cv > 0.0) {		# Add ageing error
      ## !! To do ###***
    }
    # Survey at Age
    surveyage[,i] <- agecomp
    # Survey at Length
    surveylen[,i] <- sample_lengths(ages1, lenbins, agecomp[], om$growth[[om$sex[i]]])
  }
  ##* use switch here or at least if else
  if(type=="index") return(survey)
  if(type=="age") return(round(surveyage,0))
  if(type=="len") return(round(surveylen,0))
}

#' Sampling - Catch-at-Age or Length
#'
#' Sampling - Catch-at-Age or Length
#' @param type sampling type
#' @param om operating model object
#' @param sampling sampling object
#' @param fleet fleet object
#' @param obs observations?
#' @param ff fleets???
#' @param y year???
#' @param ss ss???
#' @param rr region???
#' @export
sample_catch	<- function(type, om, sampling, fleet, obs, ff, y, ss, rr) {
  ## Sampling of catch-at-age
  ages1 <- om$ages + om$growth_props[ss]
  lenbins <- sampling$len_classes
  ageing_cv	<- sampling$ageing_error
  # Storage
  catchage <- obs[[ff]]$age_sample_n[,y,,ss,rr] * 0
  catchlen <- obs[[ff]]$len_sample_n[,y,,ss,rr] * 0
  for (i in 1:length(om$sex)) {
    catchN <- sum(fleet[[ff]]$landings_n[,y,i,ss,rr])
    if (catchN == 0) {
      ###*** should something happen here?
      #print(paste("No catch to sample for age and length (landings_n = 0) for sex",i,"in fishery",ff,"in year",y))
    } else {
      agecomp <- fleet[[ff]]$landings_n[,y,i,ss,rr] 			# Numbers at age caught
      ## Catch at Age
      if(type=="age") {
        sampleN <- as.numeric(obs[[ff]]$age_sample_N[,y,i,ss,rr])
        if (sampleN > catchN) {
          warning(paste("Catch-at-age sample N is larger than total catch N (landings_n)",ff,"in year",y,"and is reduced to catch N"))
          sampleN <- catchN
        }
        catchage[,i] <- sample_ages(agecomp, sampleN, ageing_cv)
      }
      ## Catch at Length
      if(type=="len") {
        sampleN <- as.numeric(obs[[ff]]$len_sample_N[,y,i,ss,rr])
        if (sampleN > catchN) {
          warning(paste("Catch-at-length sample N is larger than total catch N (landings_n)",ff,"in year",y,"and is reduced to catch N"))
          sampleN <- catchN
        }
        sampleage <- sample_ages(agecomp, sampleN, ageing_cv)	# Sampling catch-at-age (use agecomp)
        catchlen[,i] <- sample_lengths(ages1, lenbins, sampleage, om$growth[[om$sex[i]]])
      }
    }
  }
  ##* use switch here or at least if else
  if(type=="age") return(catchage)
  if(type=="len") return(catchlen)
}

#' Sampling - Catch rates
#'
#' Sampling - Catch rates
#' @param pop population object
#' @param sampling sampling object
#' @param fleet fleet object
#' @param ff fleets???
#' @param y year???
#' @param ss ss???
#' @param rr region???
#' @export
sample_cpue	<- function(pop, sampling, fleet, ff, y, ss, rr) {
  ## Catch rates - relative index of abundance = proportional to vulnerable biomass (by sex), arbitrarily scaled by 1000
  cpue <- apply(pop$n[,y,,ss,rr] * fleet[[ff]]$landings_sel[,y,,ss,rr] * pop$wt[,y,,ss,rr],c(2),sum) * sampling$cpue_q
  # toteff 	<- sum(fleet[[ff]]$effort[1,y,,ss,rr])	# don't use effort here
  ##* use switch here or at least if else
  if (sampling$cpue_var_real_units == "normal") cpue <- cpue * (rnorm(1,mean=1,sd=sampling$cpue_var_real))
  if (sampling$cpue_var_real_units == "lognormal") cpue <- cpue * exp(rnorm(1,mean=0,sd=sampling$cpue_var_real))
  return(cpue)
}

#' Age sample
#'
#' Age sample
#' @param agecomp	Numbers at age in the total catch
#' @param sampleN	Sample sizes for age sample
#' @param ageing_cv Ageing error
#' @export
sample_ages <- function (agecomp, sampleN, ageing_cv = 0.0) {
  ## Multinomial age samples
  nages <- length(agecomp)
  prop.age <- as.vector(agecomp/sum(agecomp))		# Proportions
  agesample <- rmultinom(n=1, size=sampleN, prob=prop.age)
  ## Add ageing error
  if(ageing_cv > 0.0) {
    ##* implement ageing error
    ## !! To do ###***
  }
  return(agesample)
}

#' Size sample
#'
#' Sample lengths
#' @param ages1 Number at age of the sample
#' @param lenbins  Length bins by year, season and region
#' @param agecomp Number at age of the sample
#' @param growth growth parameters
#' @export
sample_lengths <- function (ages1, lenbins, agecomp, growth) {
  # Mean length
  alk <- calc_VBlen(ages1, agecomp, growth)
  # Calculate normal distribution of length by length classes & sum up
  res <- vector(mode="numeric", length=length(lenbins))
  for (aa in 1:length(agecomp)) {
    # SD = Mean * CV, cumulative distribution
    res1 <- pnorm(lenbins, mean=alk[aa], sd=alk[aa]*growth[4])
    res11 <- res1[1]
    res1[1:(length(res1)-1)] <- (res1[2:length(res1)] - res1[1:(length(res1)-1)])
    res1[1] <- res1[1] + res11
    res1[length(res1)] <- 0
    res <- res + res1 * agecomp[aa]
  }
  #return(round(res,0))
  return(res)
}
