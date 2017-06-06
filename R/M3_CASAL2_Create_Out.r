#### Create output.csl file for CASAL Assessment

#' Create CASAL output.csl file
#'
#' Create CASAL output.csl file
#' @inheritParams casal
#' @export
create_casal_file_out <- function(params, casal_path, skel_csl, csl) {
	# params: 	 Object with input parameters for CASAL assessment
	# casal_path <- ""
	# skel_csl  <- "casal_output_skel.csl"
	# csl 		<- "casal_output.csl"
	# params    <- datass
  # extract casal parameters from skeleton file
	casalout <- casal::extract.csl.file(paste(casal_path, skel_csl, sep=""))
	#### 0. Prepare csl file (for sequence, see Manual Chapter 14.1.2 The input parameter files)
	# Remove some list elements (that relate to individual fisheries or tagging data)
	rm_name	<- which(substr(names(casalout),1,9) %in% c("numbers_a","abundance","n_project","quantitie"))
	casalout[rm_name]	<- NULL
	casalout[["print"]] <- NULL
	## 1. Output file: Fill in data
	casalout[["print"]] <- params$output[["print"]]
	casalout[["print"]]$command	<- "print"
	casalout[["print"]]$value	<- character(0)
	casalout[["quantities"]] <- params$output[["quantities"]]
	casalout[["quantities"]]$command <- "quantities"
	casalout[["quantities"]]$value <- character(0)
	#casalout[["abundance[1]"]]	<- params$output[["abundance[1]"]]
	for (ee in 1:length(params$regions)) {
		Natage	<- paste("numbers_at[Numbers_at_age_",params$regions[ee],"]",sep="")
		casalout[[Natage]] <- params$output[[Natage]]
	}
	casalout[["n_projections"]]$command <- "n_projections"
	casalout[["n_projections"]]$value <- params$output[["n_projections"]]
	# write the casal output file
	casal::write.csl.file(casalout,paste(casal_path,csl,sep=""))
	## Note: nothing is returned
}
