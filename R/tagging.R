## functions relating to the tagged population and tagging (moved from
## M2_Functions.R)

#' Tag binomial
#'
#' tag_binom applies to every level of dimension (by units,area,yeartag,areatag).
#' @param dat1 Data
#' @param prob1 Probabilities
#' @param dims dimensions
#' @export
tag_binom <- function(dat1, prob1, dims) {
  dat <- as.vector(round(dat1,0))
  # Expand by extra dimensions, e.g. Nyear*Narea
  probs <- rep(as.vector(prob1),dims)
  # Apply binomial distribution to all data points
  datnew <- rbinom(n = length(dat), size = dat, prob = probs)
  # Convert back
  dat1[] <- datnew
  return(dat1)
}

#' Apply natural Mortality to iTags
#'
#' Apply prob (natM) to array dat (iTags)
#' @param dat Object with individual tags
#' @param prob Object with probabilities (e.g. natural moratlity) by age, sex and area
#' @export
iTags_M <- function(dat, prob) {
  # Index of tagged fish that are alive
  liveT <- which(dat[,"Dead"]==0)		# Saveguard to exclude dead fish
  # Add M to dat2
  dat2 <- dat[liveT,]
  if(min(prob) == max(prob)) {
    # If Prob is constant across all ages, sex and areas: use first value
    dat2$prob <- min(prob)
  } else {						# if Prob is not constant: assign age, sex and area-specific value to dat2
    m_age <- rep(rep(dimnames(prob)$age,length(dimnames(prob)$sex)),length(dimnames(prob)$area))
    m_sex <- rep(rep(dimnames(prob)$sex,each=length(dimnames(prob)$age)),length(dimnames(prob)$area))
    m_area <- rep(dimnames(prob)$area,each=length(dimnames(prob)$age) * length(dimnames(prob)$sex))
    prob2 <- cbind(prob=as.vector(prob),id=paste(m_sex,"_",m_age,"_",m_area,sep=""))		# ID For matching
    dat2$id <- paste(dat2[,"Sex"],"_",dat2[,"LastAge"],"_",dat2[,"LastArea"],sep="")		# ID For matching
    dat2 <- merge(dat2,prob2,by="id",all.x=T,all.y=F)
  }
  # Apply random Binomial distribution for M to all live fish
  dead <- rbinom(n = length(liveT), size = 1, prob = as.numeric(as.character(dat2$prob)))
  # Assign dead fish to dat
  dat[liveT[dead %in% 0],"M"] <- 1
  dat[liveT[dead %in% 0],"Dead"] <- 1
  return(dat)
}

#' Release of individual tags
#'
#' Create to array dat (iTags)
#' @param dat Object with individual tags
#' @param newtags Vector of released tag numbers by age (prior to tag-release mortality & tag loss)
#' @param y Year
#' @param i Sex
#' @param ss Season
#' @param rr Region
#' @param om Object with parameters for operating model
#' @param sampling Object with parameters for sampling
#' @export
iTags_release <- function(dat, newtags, y, i, ss, rr, om, sampling) {
  ## Store individual tagged fish
  tages <- NULL
  for (x in 1:length(newtags)) tages	<- c(tages, rep(om$ages[x], newtags[x]))
  ttags <- array(0,dim=c(length(tages),dim(dat)[2]))
  colnames(ttags) <- colnames(dat)
  ttags[,"Sex"] <- om$sex[i]
  ttags[,"RelY"] <- ttags[,"LastY"] <- om$years[y]
  ttags[,"RelAge"] <- ttags[,"LastAge"] <- tages
  ttags[,"RelSeas"] <- ttags[,"LastSeas"] <- om$seasons[ss]
  ttags[,"RelArea"] <- ttags[,"LastArea"] <- om$regions[rr]
  ## Account for release M and tag loss
  dead <- rbinom(n=length(tages), size=1, prob=(1-sampling$tag_mort))
  ttags[dead %in% 0,"RelM"] <- 1
  ttags[dead %in% 0,"Dead"] <- 1
  ## Combine with dat and return
  dat <- rbind(dat, as.data.frame(ttags))
  dat[, -c(which(colnames(dat) %in% "Sex"))] <- sapply(dat[, -c(which(colnames(dat) %in% "Sex"))], as.numeric)
  return(dat)
}

#' Recapture of individual tags
#'
#' Apply prob (recapture probability) to array dat (iTags)
#' @param dat Object with individual tags
#' @param prob Recapture probability by age and sex
#' @param y Year
#' @param ss Season
#' @param rr Region
#' @param om Object with parameters for operating model
#' @export
iTags_recapture <- function(dat, prob, y, ss, rr, om) {
  # Index of tagged fish that are alive and currently in area rr
  liveT <- which(dat[,"Dead"]==0 & dat[,"LastY"]==om$years[y] &
                   dat[,"LastSeas"]==om$seasons[ss] & dat[,"LastArea"]==om$regions[rr])
  # Add age and sex-specific probs to dat2 (not exactly correct where LastAreas !=rr, but they are not used anyway)
  dat2 <- dat[liveT,]
  m_age <- rep(dimnames(prob)$age,length(dimnames(prob)$sex))
  m_sex <- rep(dimnames(prob)$sex,each=length(dimnames(prob)$age))
  prob2 <- cbind(prob=as.vector(prob),id=paste(m_sex,"_",m_age,sep=""))		# ID For matching
  dat2$id <- paste(dat2[,"Sex"],"_",dat2[,"LastAge"],sep="")					# ID For matching
  dat2 <- merge(dat2,prob2,by="id",all.x=T,all.y=F)
  # Apply random Binomial distribution for M to all live fish
  dead <- rbinom(n = length(liveT), size = 1, prob = (1-as.numeric(as.character(dat2$prob))))
  # Assign dead fish to dat
  dat[liveT[dead %in% 0],"F"] <- 1
  dat[liveT[dead %in% 0],"Dead"] <- 1
  return(dat)
}
