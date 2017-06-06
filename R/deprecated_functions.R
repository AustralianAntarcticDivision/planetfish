## Home for functions that will be removed in forthcoming versions
##
## NOTE: Add a deprecation warning to functions that are moved here

#' Initial age structure at start of first year
#'
#' Initial age structure at start of first year
#'
#' Check the definition of sigma wrt to the lognormal dist
#' @param mu mean recruitment
#' @param sigma sd/mu
#' @param age vector of (min_age, max_age)
#' @param M vector of natural mortality by age class
#' @export
init_age_comp <- function(mu, sigma, age, M) {
  warning("The function 'init_age_comp' has been replaced by 'initial_ages' and
    will be removed in version 0.7")
  n_ages <- age[2] - age[1] + 1
  # Vector with initital Recs (AveRec*sigmaRec)
  init1 <- mu*exp(rnorm(n_ages,sd=(sigma)))
  # Lower boundary: at least 100
  init1 <- ifelse(init1<(mu/100),(mu/100),init1)
  # Account for natM
  init2 <- init1[(age[1]+1):(n_ages)]*exp(-cumsum(M[age[1]:(n_ages-1)]))
  #init2[n_ages-1] <- mu*exp(-sum(M[age[1]:(n_ages-2)]))				# Assume average rec for plusgroup
  # NatM for last age
  lastM <- exp(-(M[n_ages-1]))
  initAge <- c(init1[1],init2[1:(n_ages-2)],init2[n_ages-1]*lastM/(1-lastM))
  # Return pop numbers by sex & area
  return(round(initAge,0))
}

#' Stock-recruitment relationship
#'
#' Stock-recruitment relationship
#' @param pin recruitment input parameters
#' @param mu mean recruitment
#' @param sigma sd/mu
#' @param SSB0 SSB0
#' @param SSB SSB
#' @param rec_series A recruitment series (vector)
#' @param h steepness
#' @export
SR <- function(pin, mu, sigma, SSB0 = 1, SSB = 0, rec_series = 0, h = 1) {
  warning("The function 'SR' has been replaced by 'stock_recruit' and will be
          removed in version 0.7")
  # lognormal
  if (pin == "lognormal") res <- mu * exp(rnorm(1, sd=sigma))
  # Sampling historical recruitment
  if (pin == "historical") res <- sample(rec_series,1)
  # Beverton-Holt
  if (pin == "BH") res <- (SSB/SSB0) / (1-((5*h-1)/(4*h))*(1-(SSB/SSB0))) * mu * (exp(rnorm(1,sd=sigma)))
  # Ricker
  if (pin == "Ricker") res <- (SSB/SSB0) * (1/(5*h))^((5/4)*((SSB/SSB0)-1))   * mu * (exp(rnorm(1,sd=sigma)))
  # User-defined
  if (pin == "defined") res <- rec_series
  # Recruitment must be above min limit
  res <- ifelse(res < 1000, 1000, res)
  # return the recruitment
  return(res)
}

#' Create a movement rule in the operating model
#'
#' Create a movement rule in the operating model
#'
#' The output of this function is written to om$move_rules in the original
#' implementation of movement. The columns Sex, Year and Season have been added
#' after ID. A zero is used to represent
#'
#'
#' @param move_matrix a matrix with a row for movements between each region and
#' columns with names ("Name", "Sex", "Year", "Season", "Origin", "Destination",
#' "Ages")
#' @param om operating model parameters
#' @export
create_move_rule <- function(move_matrix, om){
  warning("The function 'create_move_rule' has been replaced and will be removed
          in version 0.6")
  ##* check the function inputs
  # "Name", "Origin", "Destination", ages
  ## create a dataframe
  obj <- as.data.frame(move_matrix)
  obj[, -c(1)] <- sapply(obj[, -c(1)], function(x) as.numeric(as.character(x)))
  ## return the movement rule
  obj
}

#' Logistic ogive
#'
#' Logistic ogive
#' @param ages vector of ages
#' @param params logistic parameters
#' @export
ogive_logistic <- function(ages, params){
  warning("The individual ogive functions have been replaced by 'ogive' and will
          be removed in version 0.7")
  # 1/(1+19^((x50-ages)/x95)):
  # x50:   value of x for which ogive = 0.5
  # x95:   value to be added to x50 to give ogive = 0.95
  res <- 1/(1+19^((params[1]-ages)/params[2]))
  return(res)
}

#' Single ramp ogive
#'
#' Single ramp ogive
#' @param ages vector of ages
#' @param params logistic parameters
#' @export
ogive_ramp <- function(ages, params){
  warning("The individual ogive functions have been replaced by 'ogive' and will
          be removed in version 0.7")
  # params[1]: x value where ramp starts to increase from 0
  # params[2]: x value where ramp peaks at 1
  up <- params[2]-params[1] +1
  res <- c(rep(0,params[1]-ages[1]),round(seq(0,1,length.out=up),4),
           rep(1,ages[length(ages)]-params[2]))
  return(res)
}

#' Double gaussian ogive
#'
#' Double gaussian ogive
#' @param ages vector of ages
#' @param params logistic parameters
#' @export
ogive_dbnormal <- function(ages, params) {
  warning("The individual ogive functions have been replaced by 'ogive' and will
          be removed in version 0.7")
  # x <= a: 2^{-[(ages-a)/sigma.left]^2}	; x > a: 2^{-[(ages-a)/sigma.right]^2}
  # params[1]: a = top
  # params[2]: sigma.left
  # params[3]: sigma.right
  res <- 2^(-(((ages-params[1])/ifelse(ages <= params[1],params[2],params[3]))^2))
  return(res)
}

#' Double plateau gaussian ogive
#'
#' Double plateau gaussian ogive
#' @param ages vector of ages
#' @param params logistic parameters
#' @export
ogive_dbplateaunormal <- function(ages, params) {
  warning("The individual ogive functions have been replaced by 'ogive' and will
          be removed in version 0.7")
  # Input parameters (see Bull et al. 2005 - CASAL, p.48):
  # x <= a1: 2^[(ages-a1)/sigma.left]^2;
  # a1 < x <= a1+a2: a.max;
  # x > a1+a2: 2^[(ages-(a1+a2))/sigma.right]^2
  # params[1]: a1 = lower.bound for plateau
  # params[2]: a2 = upper.bound for plateau
  # params[3]: sigma.left
  # params[4]: sigma.right
  # a.max set = 1
  a1 <- params[1]
  a2 <- params[2]
  res <- ifelse(ages <= a1, 2^(-(((ages-a1)/params[3])^2)),
                ifelse(ages > (a1+a2), 2^(-(((ages-(a1+a2))/params[4])^2)),1))
  return(res)
}

#' Double ramp ogive
#'
#' Double ramp ogive
#' @param ages vector of ages
#' @param params logistic parameters
#' @export
ogive_dbramp <- function(ages, params){
  warning("The individual ogive functions have been replaced by 'ogive' and will
          be removed in version 0.7")
  # ages: vector for which logistic values are calculated
  # params[1]: low 0 -  x value where 1st ramp starts to increase from 0
  # params[2]: low 1 -  x value where 1st ramp peaks at 1
  # params[3]: high 1 - x value where 2nd ramp starts to decrease from 1
  # params[4]: high 0 - x value where 2nd ramp reaches 0
  up   <- params[2]-params[1] + 1
  down <- params[4]-params[3] + 1
  res <- c(rep(0,params[1]-ages[1]),round(seq(0,1,length.out=up),4),rep(1,(params[3]-params[2]-1)),
           round(seq(1,0,length.out=down),4),rep(0,ages[length(ages)]-params[4]))
  return(res)
}

#' Provide input data for ogive
#'
#' Provide input data for ogive
#' @param ages vector of ages
#' @param params ogive parameters
#' @export
ogive_provide <- function(ages,params){
  warning("The individual ogive functions have been replaced by 'ogive' and will
          be removed in version 0.7")
  params
}

#' Movement function
#'
#' The original movement function converted from a matrix to a dataframe
#'
#' #' Note this function should only be used when movement does not vary with Sex,
#' Year or Season
#'
#' * use a better selection method than dimension of array
#' it should be a switch with a check
#' @param dat Population or tag numbers
#' @param om operating model object
#' @export
move_fish_dframe_old <- function (dat, om)  {
  ## checks
  if(any(om$move_rules[["Sex"]] != 0))
    stop("The movement rule specifies sex specific movement which is not
         implemented for dframe movement")
  if(any(om$move_rules[["Year"]] != 0))
    stop("The movement rule specifies year specific movement which is not
         implemented for dframe movement")
  if(any(om$move_rules[["Season"]] != 0))
    stop("The movement rule specifies season specific movement which is not
         implemented for dframe movement")
  # dat 	<- pop$n[,y,,ss,,drop=FALSE]
  # dat	<- tag$tags[,y,,ss,,,,drop=FALSE]
  # Get movement rules
  mrule <- om$move_rules
  # First create list with numbers of fish that move out of an area
  move_n <- list()
  ## calculate the number of fish moving
  for (i in 1:nrow(mrule)) {
    if(length(dim(dat)) == 5) 		# pop$n
      move_n[[i]] <- dat[,,,,mrule$Origin[i]] * as.numeric(mrule[mrule$Origin[i],colnames(mrule)%in%om$names_ages])
    if(length(dim(dat)) == 7) 		# tag$tags
      move_n[[i]] <- dat[,,,,mrule$Origin[i],,] * as.numeric(mrule[mrule$Origin[i],colnames(mrule)%in%om$names_ages])
  }
  ## move the fish
  for (i in 1:nrow(mrule)) {
    if(length(dim(dat)) == 5) {		# pop$n
      dat[,,,,mrule$Origin[i]] <- dat[,,,,mrule$Origin[i]] - move_n[[i]]
      dat[,,,,mrule$Destination[i]] <- dat[,,,,mrule$Destination[i]] + move_n[[i]]
    }
    if(length(dim(dat)) == 7) {	# tag$tags
      dat[,,,,mrule$Origin[i],,] <- dat[,,,,mrule$Origin[i],,] - move_n[[i]]
      dat[,,,,mrule$Destination[i],,] <- dat[,,,,mrule$Destination[i],,] + move_n[[i]]
    }
  }
  ## return the population array
  return(dat)
}

