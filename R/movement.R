## Movement functions live here
## Movement could do with its own help file
## There is also a movement markdown document, make into a vignette

#' Initial Movement
#'
#' Apply movement to create initial population numbers by age, sex and area
#'
#' Note this function should only be used when movement does not vary with Sex,
#' Year or Season
#' @param n Initial population numbers
#' @param om operating model object
#' @param equilib_yrs number of years to move the population prior to the first
#' model timestep (default=100).
#' @export
move_init <- function (n, om, equilib_yrs = 100)  {
  ## checks
  if(any(om$move_rules[["Sex"]] != 0))
    stop("The movement rule specifies sex specific movement which is not
         implemented for initial movement")
  if(any(om$move_rules[["Year"]] != 0))
    stop("The movement rule specifies year specific movement which is not
         implemented for initial movement")
  if(any(om$move_rules[["Season"]] != 0))
    stop("The movement rule specifies season specific movement which is not
         implemented for initial movement")
  ##* I think this may need to be over sex
  #n <- mod$init_age # remember to disable
  mrule <- om$move_rules
  for (x in 1:equilib_yrs) {
    # Ageing
    n[om$n_ages,,,,] <- (n[om$n_ages,,,,] + n[om$n_ages-1,,,,]) * exp(-om$natM[om$n_ages])
    n[2:(om$n_ages-1),,,,] <- n[1:(om$n_ages-2),,,,] * exp(-om$natM[1:(om$n_ages-2)])
    recs <- om$rec_mu * exp(rnorm(1,sd=(om$rec_sigma)))
    ## assign recruitment to regions
    for (rr in 1:om$n_regions){
      n[1,1,,1,rr] <- recs * om$rec_area[rr] * om$rec_sex
    }
    # Movement
    # First create list with numbers of fish that move out of an area
    move_n	<- list()
    for (i in 1:nrow(mrule)){
      move_n[[i]] <- n[,,,,mrule$Origin[i]] *
        as.numeric(mrule[mrule$Origin[i], colnames(mrule) %in% om$names_ages])
    }
    # Then move fish
    for (i in 1:nrow(mrule)) {
      n[,,,,mrule$Origin[i]] <- n[,,,,mrule$Origin[i]] - move_n[[i]]
      n[,,,,mrule$Destination[i]] <- n[,,,,mrule$Destination[i]] + move_n[[i]]
    }
  }
  return(n)
}

#' Movement function
#'
#' The original movement function converted from a matrix to a dataframe
#'
#' Note this function should only be used when movement does not vary with Sex,
#' Year or Season
#' @param n Population or tag numbers
#' @param move_rule movement dataframe specifying movement by region, sex, year,
#' season and age classes
#' @param age_names vector of age class names
#' @param type switch to specifiy movement of population ("pop") or tags ("tag")
#' @export
move_fish_dframe <- function (n, move_rule, age_names, type="pop")  {
  ## checks
  if(any(move_rule[["Sex"]] != 0))
    stop("The movement rule specifies sex specific movement which is not
         implemented for dframe movement")
  if(any(move_rule[["Year"]] != 0))
    stop("The movement rule specifies year specific movement which is not
         implemented for dframe movement")
  if(any(move_rule[["Season"]] != 0))
    stop("The movement rule specifies season specific movement which is not
         implemented for dframe movement")
  # First create list with numbers of fish that move out of an area
  move_n <- list()
  # move population or tags
  switch(type,
         ## population movement
         pop = {if(length(dim(n)) != 5) stop("population array has the wrong
                                               dimensions")
           ## calculate the numbers moving by age
           for (i in 1:nrow(move_rule)) {
             move_n[[i]] <- n[,,,,move_rule$Origin[i]] *
               as.numeric(move_rule[move_rule$Origin[i],colnames(move_rule) %in% age_names])
           }
           ## move the population
           for (i in 1:nrow(move_rule)) {
             n[,,,,move_rule$Origin[i]] <- n[,,,,move_rule$Origin[i]] - move_n[[i]]
             n[,,,,move_rule$Destination[i]] <- n[,,,,move_rule$Destination[i]] +
               move_n[[i]]
           }
         },
         ## tag movement
         tag = {if(length(dim(n)) != 7) stop("tag array has the wrong
                                               dimensions")
           ## calculate the numbers moving by age
           for (i in 1:nrow(move_rule)) {
             move_n[[i]] <- n[,,,,move_rule$Origin[i],,] *
               as.numeric(move_rule[move_rule$Origin[i],colnames(move_rule) %in% age_names])
           }
           ## move the tagged fish
           for (i in 1:nrow(move_rule)) {
             n[,,,,move_rule$Origin[i],,] <- n[,,,,move_rule$Origin[i],,] - move_n[[i]]
             n[,,,,move_rule$Destination[i],,] <- n[,,,,move_rule$Destination[i],,] +
               move_n[[i]]
           }
         })
  ## return the population array
  n
}

#' Movement function for individual tags
#'
#' Movement of individual tags, based on move_rules
#'
#' This implementation of individual based tag movement is designed for use when
#' a large number of regions are specified and array based movement becomes slow
#' @param dat Object with individual tags
#' @param om operating model object
#' @export
move_iTags	<- function(dat, om) {
  # dat = tag$iTags
  # om: Object with om parameters
  # xtabs(~Dead + LastArea, data=dat)

  # Get movement rules / probabilities
  mrule <- om$move_rules		# Get movement rules
  for (i in 1:nrow(mrule)) {
    # Index of tagged fish that are alive and currently in area
    liveT <- which(dat[,"Dead"]==0 & dat[,"LastArea"] == mrule$Origin[i])
    mrule1 <- t(rbind(prob=mrule[i,colnames(mrule)%in%om$names_ages],Age=om$ages))
    # Select fish and add movement prob
    dat2 <- merge(dat[liveT,],mrule1,by.x="LastAge",by.y="Age",all.x=T,all.y=F)
    # xtabs(~prob + LastAge, data=dat2)
    mm <- rbinom(n = length(liveT), size = 1, prob = as.numeric(as.character(dat2$prob)))
    dat[liveT[mm %in% 1],"LastArea"]	<- mrule$Destination[i]
  }
  return(dat)
}


#' Expand a movement matrix
#'
#' Expand a movement matrix.
#'
#' This function replaces 'create_move_rule' for array based movement. While the
#' current implementation is slow it is only used once in every iteration.
#' @param move_matrix a matrix with a row for movements between each region and
#' columns with names ("Name", "Sex", "Year", "Season", "Origin", "Destination",
#' "Ages")
#' @param om operating model parameters
#' @export
expand_move_matrix <- function(move_matrix, om){
  ## define an empty dataframe for the results
  res <- data.frame("Origin" = numeric(0),"Destination" = numeric(0),
                    "Sex" =  numeric(0), "Year" =numeric(0), "Season" = numeric(0),
                    setNames(replicate(om$n_ages,numeric(0), simplify = F), om$names_ages))
  ## loop over regions, sexes, years, seasons
  for(orig_Region in 1:om$n_regions){
    for(dest_Region in 1:om$n_regions){
      for(iSex in 1:om$n_sexes){
        for(iYear in 1:om$n_years){
          for(iSeason in 1:om$n_seasons){
            ## select the ith movement rule
            i_move <- move_matrix[move_matrix[["Origin"]] == orig_Region &
                            move_matrix[["Destination"]] == dest_Region &
                            move_matrix[["Sex"]] %in% c(iSex, 0) &
                            move_matrix[["Year"]] %in% c(iYear, 0) &
                            move_matrix[["Season"]] %in% c(iSeason, 0),
                          colnames(move_matrix) %in% om$names_ages]
            ## no movement fill with zeroes
            if(nrow(i_move) == 0){
              i_row <- data.frame("Origin" = orig_Region, "Destination" = dest_Region,
                                  "Sex" = iSex, "Year" = iYear, "Season" = iSeason,
                                  setNames(replicate(om$n_ages, 0, simplify = F),
                                           om$names_ages))
              res <- rbind(res, i_row)
              ## otherwise specify the movement by age class and other covariates
            }else if(nrow(i_move) == 1){
              i_row <- data.frame("Origin" = orig_Region, "Destination" = dest_Region,
                                  "Sex" = iSex, "Year" = iYear, "Season" = iSeason,
                                  setNames(i_move, om$names_ages))
              ## append to the res
              res <- rbind(res, i_row)
            }else stop("More than one movement rule is specified for an area, sex,
                       year and season")
          }
        }
      }
    }
  }
  ## return the results
  res
}
