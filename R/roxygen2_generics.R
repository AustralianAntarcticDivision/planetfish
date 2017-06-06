# some generic documentation

#' Operating model parameters
#'
#' @param para list of model parameters.
#' @param res list to store model results.
#' @param tag tagging data
#' @name om
NULL
#> NULL

#' CASAL assessment parameters
#'
#' @param params assessment parameters (could replace with datass)
#' @param datass assessment data
#' @param tag tagging data
#' @param casal_path path to casal.exe
#' @param mpd_dat mpd data
#' @param mvnsamples_dat mvn samples
#' @param inputprefix CASAL input file prefix
#' @param output_log output log file
#' @param linux run on Windows (default 0) or Linux (1)
#' @param intern CASAL is run internally (default = TRUE) with no output sent
#' to the R console or it is run with output sent to R console (FALSE, use for
#' debugging)
#' @param skel_csl skeleton csl file (or location of file)
#' @param csl csl file (name perhaps?)
#' @param pop_csl CASAL population csl file (name, location???)
#' @param TAC Initial catch level
#' @param use_specific_ratios If 1, then use TAC ratios provided by 'specific_ratios'
#' @param specific_ratios Specific TAC ratios (if use_specific_ratios=1)
#' @param Yr_current current year
#' @param ssb_level_deplet depletion level
#' @param ssb_level_target target level
#' @param ref_prob_deplet reference depletion probability
#' @param ref_prob_target reference target probability
#' @param tol tolerance
#' @param set_seed switch for random number seed
#' @param rand_seed random seed
#' @name casal
NULL
#> NULL





