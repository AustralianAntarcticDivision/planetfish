#' Ogive function
#'
#' Provide input data for ogive
#'
#' Add extended details to the help file
#'
#' logistic =
#'
#' ramp =
#'
#' dbnormal =
#'
#' dbnormalplateau = Input parameters (see Bull et al. 2005 - CASAL, p.48):
#' x <= a1: 2^[(ages-a1)/sigma.left]^2;
#' a1 < x <= a1+a2: a.max;
#' x > a1+a2: 2^[(ages-(a1+a2))/sigma.right]^2
#' params[1]: a1 = lower.bound for plateau
#' params[2]: a2 = upper.bound for plateau
#' params[3]: sigma_left
#' params[4]: sigma_right
#' a.max set = 1
#'
#' dbramp = # params[1]: low 0 -  x value where 1st ramp starts to increase from 0
#' params[2]: low 1 -  x value where 1st ramp peaks at 1
#' params[3]: high 1 - x value where 2nd ramp starts to decrease from 1
#' params[4]: high 0 - x value where 2nd ramp reaches 0
#' @param type ogive type, either "logistic", "ramp", "dbnormal",
#' "dbplateaunormal", "dbramp", "provide"
#' @param ages vector of ages
#' @param params ogive parameters (vary with type)
#' @export
ogive <- function(type, ages, params){
  ## checks
  ## define valid ogives, remember to add new ogives below
  valid_ogives <- c("logistic", "ramp", "dbnormal", "dbplateaunormal", "dbramp",
                    "provide")
  ## do some checks on res
  if(!type %in% valid_ogives)
    stop(paste0("supplied ogive type must be one of ", valid_ogives))
  ## check the input parameters
  # switch(type,
  #        # !any(true_vars %in% names(res))
  #        ##* below isn't quite right
  #        logistic = if(!any(params %in% c("x50", "x95"))),
  ## calculate proportions by age based on specified type
  switch(type,
         logistic = res <- 1/(1+19^((params[["x50"]]-ages)/params[["x95"]])),
         ramp = res <- c(rep(0,params[["start"]]-ages[1]),
                         round(seq(0,1,length.out=(params[["peak"]] -
                                                     params[["start"]] + 1)),4),
                         rep(1,ages[length(ages)]-params[["peak"]])),
         dbnormal = res <- 2^(-(((ages-params[["top"]])/ifelse(ages <= params[["top"]],params[["sigma_left"]],params[["sigma_right"]]))^2)),
         dbplateaunormal = res <- ifelse(ages <= params[["a1"]],
                                         2^(-(((ages-params[["a1"]])/params[["sigma_left"]])^2)),
                                         ifelse(ages > (params[["a1"]]+params[["a2"]]),
                                                2^(-(((ages-(params[["a1"]]+params[["a2"]]))/params[["sigma_right"]])^2)),1)),
         dbramp = {
           # params[1]: low0 -  x value where 1st ramp starts to increase from 0
           # params[2]: low1 -  x value where 1st ramp peaks at 1
           # params[3]: high1 - x value where 2nd ramp starts to decrease from 1
           # params[4]: high0 - x value where 2nd ramp reaches 0
           up <- params[["low1"]]-params[["low0"]] + 1
           down <- params[4]-params[3] + 1
           res <- c(rep(0,params[["low0"]]-ages[1]),round(seq(0,1,length.out=up),4),rep(1,(params[3]-params[["low1"]]-1)),
                    round(seq(1,0,length.out=down),4),rep(0,ages[length(ages)]-params[4]))},
         provide = res <- params)
  ## return the ogive
  res
}
