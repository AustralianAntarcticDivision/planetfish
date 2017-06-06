## process outputs from simulations

#' Extract simulation data and calculate relative error
#'
#' Extract simulation data and calculate relative error
#' @param sim_data list of outputs from simulations (currently csv files)
#' @param var_names list of operating model and assessment model parameter names
#' to calculate relative error from
#' @export
rel_err <- function(sim_data, var_names){
  ## create an empty list with the variable names
  res <- vector("list", length(var_names[["am"]]))
  names(res) <- var_names[["am"]]
  ## loop over the variables to estimate relative error for each scenario
  for(i in 1:length(var_names$om)){
    ## object to store the estimates for each parameter
    scenarios  <- vector("list", length(names(sim_data)))
    names(scenarios) <- names(sim_data)
    for(j in 1:length(sim_data)){
      ## calculate rel err = est - true / true
      scenarios[[j]] <- (sim_data[[j]][[var_names$am[i]]] - sim_data[[j]][[var_names$om[i]]])/
        sim_data[[j]][[var_names$om[i]]]
    }
    ## add to the results list
    res[[i]] <- scenarios
  }
  ## return the list of relative error
  res
}
