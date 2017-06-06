#' planetfish: A fisheries operating model integreted with the CASAL
#' stock assessment model
#'
#' DRAFT
#'
#' Based on planetfish ...
#'
#' Simulation model
#'
#' The simulation model ...
#'
#' Run Annual Operating model loops with/without assessment
#'
#' Time sequence - add a link to the function that runs the model
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
#' 6. Sampling (Age samples, length samples
#'
#' 7. Tag release and recaptures
#'#'
#' Assessment model
#'
#' The CASAL integrated assessment is used ...
#'
#' CASAL sequence within a time step:
#'
#' Ageing (in an age-based model)
#'
#' Recruitment
#'
#' Maturation (if maturity is a character in the partition)
#'
#' Migration (if the model includes more than one area).
#'
#' Growth (in a size-based model).
#'
#' Natural and fishing mortality.
#'
#' Disease mortality.
#'
#' Tag release events.
#'
#' Tag shedding rate.
#'
#' Semelparous mortality.
#'
#'
#' Sampling matrix: Define data types and years by fisheries
#' for the assessment
#'
#' Data collection system
#'
#' Determine years of sampling (not season specific - assumes
#'  that all seasons with a fishery are sampled!):
#'
#' 1. catchage:  Age composition samples
#'
#' 2. catchlen:  Length composition samples
#'
#' 3. surveyage: Age survey
#'
#' 4. catchlen:  Length survey
#'
#' 5. index:     Index of abundance
#'
#' 6. tagging:   Years with tag-releases
#'
#' @docType package
#' @name planetfish
#' @importFrom stats median pnorm rbinom rmultinom rnorm runif setNames xtabs
#' @importFrom utils read.table write.table
NULL
#> NULL
