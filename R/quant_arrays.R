## create arrays for storage
## be nice to have a generic name for the lists that store the information

#' Create quant array
#'
#' Quant array
#'
#' The arrays are used to store the model objects
#' @param om list object storing the ages, lengths etc
#' @param sampling can we remove this argument
#' @param type switch
#' @param fill value to fill array with (default=0)
#' @export
create_array <- function(om, sampling, type, fill=0){
  ## add checks, might need different checks for types

  ## use switch to create the different options, names below
  switch(type,
         simple = {
           obj <- array(fill, dim=c(1,om$n_years,length(om$sex),om$n_seasons,om$n_regions),
                                   dimnames=list(quant=1,year=om$years,sex=om$sex,season=om$seasons,area=om$regions))
         },
         age = {
           obj <- array(fill, dim=c(om$n_ages,om$n_years,length(om$sex),om$n_seasons,om$n_regions),
                                 dimnames=list(age=om$ages,year=om$years,sex=om$sex,season=om$seasons,area=om$regions))
         },
         length = {
           obj <- array(fill, dim=c(sampling$n_lengths,om$n_years,length(om$sex),om$n_seasons,om$n_regions),
                                                  dimnames=list(len=sampling$len_classes,year=om$years,sex=om$sex,season=om$seasons,area=om$regions))
         },
         tag_age = {
           obj <- array(fill, dim=c(om$n_ages,om$n_years,length(om$sex),om$n_seasons,om$n_regions,
                                             om$n_years,om$n_regions),
                                    dimnames=list(age=om$ages,year=om$years,sex=om$sex,season=om$seasons,area=om$regions,
                                                  tagyear=om$years,tagarea=om$regions))
         },
         tag_length = {
           obj <- array(fill, dim=c(sampling$n_lengths,om$n_years,length(om$sex),om$n_seasons,om$n_regions,
                                                             om$n_years,om$n_regions),
                                                    dimnames=list(len=sampling$len_classes,year=om$years,sex=om$sex,season=om$seasons,area=om$regions,
                                                                  tagyear=om$years,tagarea=om$regions))
         },
         fishery = {
           obj <- array(data=fill, dim=c(om$n_years,om$n_fisheries,om$n_seasons,om$n_regions),
                        dimnames=list("Year"=om$years, "Fishery"=om$fishery, "Season"=om$seasons, "Region"=om$regions))
         },
         move = {
             obj <- array(fill, dim=c(om$n_ages,om$n_years,length(om$sex),
                                   om$n_seasons,om$n_regions, om$n_regions),
                          dimnames=list(age=om$ages,year=om$years,sex=om$sex,
                                        season=om$seasons,origin_area=om$regions,
                                        destination_area=om$regions))

           })
  ## consider adding S3 class
  ## return the array
  obj
}
