#' Assessment sampling parameters
#'
#' Assessment sampling parameters
#'
#' Determine years of sampling (not season specific - assumes that all seasons
#' with a fishery are sampled!):
#'
#' 1. catchage:  Age composition samples
#' 2. catchlen:  Length composition samples
#' 3. surveyage: Age survey
#' 4. catchlen:  Length survey
#' 5. index:     Index of abundance
#' 5. tagging:   Years with tag-releases
#'
#' sample_years is set up for all years of the OM. No need to define season & region
#' (they are fishery-specific)
#'
#' Each data type probably needs to be a list and looped over them
#' @param fleets list of the fleets (sub-fisheries)
#' @param years vector of assessment model years
#' @param types vector of sample data types
#' @param ycurr current year
#' @param age_y vector of age years??? (not sure this is used)
#' @param surv_y survey years???
#' @param catchlen_yrs years of catch at length
#' @param catchage_yrs years of catch at age
#' @param survey_yrs years of surveys
#' @param surveyage_yrs years of survey age data
#' @param surveylen_yrs years of survey length data
#' @param cpue_yrs years of CPUE data
#' @param tagging_yrs years of tagging data
#' @param cpue_cvpar Variance estimates used in assessment (default = 0.234)
#' @param cpue_cv_dist Error distribution of estimated abundance index in units
#' used in assessment either "lognormal" (default) or "normal"
#' @param cpue_cv_process_error is process error applied to CPUE (default = 0 = no,
#' 1 = yes)
#' @param surveyindexvar survey index variance (default = 0.3)
#' @param surveyagevar survey age variance (default = 0.3)
#' @param surveylenvar survey length variance (default = 0.0)
#' @param survey_dist Error distribution of survey-at-age and survey-at-length
#' data (default = "lognormal")
#' @param survey_cv_process_error	If 0, no process error is applied (= no down-weighting)
## Catch-at-length and catch-at-age
#' @param catch_at_dist	<- "multinomial" 	# Error distribution of catch-at-age and catch-at-length data
#' @param catch_at_r Robustification parameter used in binomial likelihood
#' (default = 1e-11, relevant when fitted prop is close to 1 or 0)
#' @param pin_diff_N_in_assessment If 1, use these N for assessment, if 0, use N
#'  from sampling (in OM) in assessment
#' @param catchage_N Sample size for catch at age, by sex (default = 100)
#' @param catchlen_N Sample size for catch at length, by sex (default = 100)
#' @param proportion_mortality For abundance, survey and tag-recapture
#' observations (default = 0.5)
#' @param ... additional parameters (not yet implemented)
#' @export
am_sampling <- function(fleets = list("LL"),
                        years = 1990:2010,
                        types = c("catchage","catchlen","survey","surveyage",
                                  "surveylen","cpue","tagging"),
                        ycurr = 2010,
                        age_y	= 2000:2010,
                        surv_y = seq(2000, 2010, 2),
                        ## I'm specifying the years these data sources are
                        ## available in the assessment model
                        catchlen_yrs = NULL,
                        catchage_yrs = 2000:2010,
                        survey_yrs = NULL,
                        surveyage_yrs = NULL,
                        surveylen_yrs = NULL,
                        cpue_yrs = NULL,
                        tagging_yrs = 2005:2009,
                        cpue_cvpar = 0.234,
                        cpue_cv_dist = "lognormal",
                        cpue_cv_process_error = 0,
                        surveyindexvar = 0.3,
                        surveyagevar = 0.3,
                        surveylenvar = 0.0,
                        survey_dist = "lognormal",
                        survey_cv_process_error	= 0,
                        catch_at_dist	= "multinomial",
                        catch_at_r = 1e-11,
                        pin_diff_N_in_assessment = 1,
                        catchage_N = list("LL" = 100),
                        catchlen_N = list("LL" = 100),
                        proportion_mortality = 0.5,
                        ...){
  ## create an array to hold the sampling
  first_yr <- years[1]
  ## might need to rearrange the years
  sample_years <- array(0, dim=c(length(types),length(years),length(fleets)),
                        dimnames=list(types=types,year=years,fishery=fleets))
  ## fill the array with the sample years (1 for data, 0 no data)
  if(!is.null(catchlen_yrs)) sample_years["catchlen", catchlen_yrs - first_yr + 1,] <- 1
  if(!is.null(catchage_yrs)) sample_years["catchage", catchage_yrs - first_yr + 1,] <- 1
  if(!is.null(survey_yrs)) sample_years["survey", survey_yrs - first_yr + 1, 1] <- 1
  if(!is.null(surveyage_yrs)) sample_years["surveyage", surveyage_yrs - first_yr + 1,]	<- 1
  if(!is.null(surveylen_yrs)) sample_years["surveylen", surveylen_yrs - first_yr + 1,] <- 1
  if(!is.null(cpue_yrs)) sample_years["cpue", cpue_yrs - first_yr + 1,] <- 1
  if(!is.null(tagging_yrs)) sample_years["tagging", tagging_yrs - first_yr + 1,]	<- 1
  ## now some type of loop
  #### Abundance index, Survey data, Catch-at-age and catch-at-length
  ##** these parameters are hard coded and should be set in the function args
  #   for (ff in 1:length(list_fishery)) {
  #     ##** not sure this is needed?
  #     fish	<- list_fishery[ff]	# Assessment fishery
  #   }
  ## construct list object of sampling
  sampling <- list("sample_years" = sample_years,
                   "types" = types,
                   "ycurr" = ycurr,
                   "age_y"	= age_y,
                   "surv_y" = surv_y,
                   "cpue_cvpar" = cpue_cvpar,
                   "cpue_cv_dist" = cpue_cv_dist,
                   "cpue_cv_process_error" = cpue_cv_process_error,
                   "surveyindexvar" = surveyindexvar,
                   "surveyagevar" = surveyagevar,
                   "surveylenvar" = surveylenvar,
                   "survey_dist" = survey_dist,
                   "survey_cv_process_error" =	survey_cv_process_error,
                   "catch_at_dist" = catch_at_dist,
                   "catch_at_r" = catch_at_r,
                   "pin_diff_N_in_assessment" = pin_diff_N_in_assessment,
                   "catchage_N" = catchage_N,
                   "catchlen_N" = catchlen_N,
                   "proportion_mortality" = proportion_mortality)
  ## return the sampling parameters
  sampling
}
