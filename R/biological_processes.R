## biological processes, growth, recruitment etc (moved from M2_Functions.R)

#' Von Bertalanffy length
#'
#' Growth may vary seasonally, therefore the whole
#' object is passed on
#'
#' NOTE: params1 should reference names not position
#' @param age1 age
#' @param object object to store results
#' @param params1 von Bert parameters
#' @export
calc_VBlen <- function(age1,object,params1) {
  Linf <- params1[1]; K <- params1[2]; t0 <- params1[3]; CV <- params1[4]  # CV not used here
  res <- object
  res[] <- Linf*(1-exp(-K*(age1-t0)))
  return(res)
}

#' Von Bertalanffy weight
#'
#' Growth may vary seasonally, therefore the whole (FLQuant-like)
#' object is passed on
#'
#' NOTE: params1 and params2 should reference names not position
#' @param age1 age
#' @param object object to store results
#' @param params1 von Bert parameters
#' @param params2 weight length parameters
#' @export
calc_VBweight <- function(age1,object,params1,params2) {
  Linf <- params1[1]; K   <- params1[2]; t0 <- params1[3]; CV <- params1[4]  # CV not used here
  WLa  <- params2[["a"]]; WLb <- params2[["b"]]
  res <- object
  res[] <- WLa*(Linf*(1-exp(-K*(age1-t0))))^WLb
  return(res)
}

#' Stock-recruitment relationship
#'
#' Stock-recruitment relationship
#' @param method recruitment input parameters
#' @param mu mean recruitment
#' @param sigma sd/mu
#' @param SSB0 SSB0
#' @param SSB SSB
#' @param rec_series A recruitment series (vector)
#' @param h steepness
#' @param year model year (if
#' @param rec_min minimum recruitment (default=1000)
#' @export
stock_recruit <- function(method, mu, sigma, SSB0 = 1, SSB = 0, rec_series = 0,
                          h = 1, year, rec_min = 1000) {
  ## Define mod$ssb0
  switch(method,
         lognormal = res <- mu * exp(rnorm(1, sd=sigma)),
         historical = res <- sample(rec_series,1),
         BH = res <- (SSB/SSB0) / (1-((5*h-1)/(4*h))*(1-(SSB/SSB0))) * mu * (exp(rnorm(1,sd=sigma))),
         Ricker = res <- (SSB/SSB0) * (1/(5*h))^((5/4)*((SSB/SSB0)-1)) * mu * (exp(rnorm(1,sd=sigma))),
         provide = res <- rec_series[year])
  # check the recruitment must be above min limit
  res <- ifelse(res < rec_min, rec_min, res)
  # return the recruitment
  return(res)
}

#' Initial age structure at start of first year
#'
#' Initial age structure at start of first year
#'
#' This function replaces 'init_age_comp' which will be removed in version 0.7.
#'
#' Check the definition of sigma wrt to the lognormal dist
#' @param method stoch - stochastic recruitment (standard), provide -
#' provide and include natural mortality, provide_nat without natural mortality
#' @param mu mean recruitment
#' @param sigma sd/mu
#' @param age vector of (min_age, max_age)
#' @param M vector of natural mortality by age class
#' @param age_series provide a series of age classes (methods 'provide' and
#' 'user')
#' @param R_min minimum recruitment (default=100)
#' @export
#' @examples
#' # The "stoch" method used for stochastic simulations
#' initial_ages(method="stoch", mu=1.5e+07, sigma=1e-01, age=c(1,30), M=rep(0.155,30))
#'
#' # Alternately the recruitment series is specified and natural mortality applied
#' initial_ages(method="provide", age=c(1,30), M=rep(0.155,30), age_series = rep(10000,30))
#'
#' # Or recruitment is specified without natural mortality
#' initial_ages(method="user", age=c(1,10), age_series = rep(1000,10))
initial_ages <- function(method, mu = NULL, sigma = NULL, age, M = NULL,
                         age_series = NULL, R_min=100) {
  ## calculate the number of age classes
  n_ages <- age[2] - age[1] + 1
  ## some checks
  if(method %in% c("stoch", "provide") & length(M) != n_ages)
    stop("Vector of natural mortality must have length equal to the number of
         age-classes")
  if(method %in% c("provide", "user") & length(age_series) != n_ages)
    stop("age series vector must have length equal to the number of age-classes")
  ## calculate the initial age structure based on the method
  switch(method,
         # Vector with initital Recs (AveRec*sigmaRec)
         stoch = {
           init_recs <- mu * exp(rnorm(n_ages, sd=(sigma)))
           # Lower boundary: at least
           init_recs <- ifelse(init_recs < (mu/R_min),(mu/R_min), init_recs)
         },
         ## this option includes natural mortality
         provide = {
           init_recs <- age_series
         },
         ## this
         ## perhaps have a return in here
         user = return(round(age_series,0))
  )
  # Account for natM
  adj_recs <- init_recs[(age[1]+1):(n_ages)]*exp(-cumsum(M[age[1]:(n_ages-1)]))
  # adj_recs[n_ages-1] <- mu*exp(-sum(M[age[1]:(n_ages-2)]))
  # Assume average rec for plusgroup
  last_M <- exp(-(M[n_ages-1])) # NatM for last age
  init_age <- c(init_recs[1], # first age class
                adj_recs[1:(n_ages-2)], # second to n-1 age class
                adj_recs[n_ages-1]*last_M/(1-last_M)) # last age class
  # Return pop numbers by sex & area
  return(round(init_age,0))
}


