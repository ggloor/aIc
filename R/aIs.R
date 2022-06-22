#' \code{aIc.singular} tests for singular data. 
#'   This is expected to be true if the transform is behaving rationally in
#'   compositional datasets and also true in the case of datasets with more 
#'   features than samples. 
#'
#' @param data can be any dataframe or matrix with samples by column
#' @param norm.method can be prop, clr, RLE, TMM, TMMwsp
#' @param zero.method is a logical. Filter data to remove features that are 0
#'   across all samples. Default is TRUE
#' @param log is a logical. log transform the RLE or TMM outputs, default=FALSE
#' @param group is a vector containing group information. Required fro RLE and 
#'   TMM based normalizations.
#'
#' @return Returns a list with a yes/no binary decision in 
#'   \code{is.singular} and the covariance matrix in \code{cov.matrix}
#'
#' @author Greg Gloor
#'
#' @export aIc.singular
#'
#' @importFrom matrixcalc is.singular.matrix
#'
#' @examples
#' library(ALDEx2)
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.singular(selex, norm.method='prop')

aIc.singular <- function(data, norm.method='prop', zero.method='remove', log=FALSE, group=NULL){
  
  # remove features with 0 counts across all samples only
  if(zero.method == 'remove'){
    data <- data[rowSums(data) > 0,]
  }

  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  
	cov.all <- as.matrix(cov(t(x.1)))
	
	is.singular <- matrixcalc::is.singular.matrix(cov.all)
	
	if(is.singular){
	  is.singular='Yes'
	}else {
	  is.singular='No'
	}
	return( list(is.singular=is.singular, cov.matrix=cov.all))
}
