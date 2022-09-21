#' \code{aIc.singular} tests for singular data. 
#'   This is expected to be true if the transform is behaving rationally in
#'   compositional datasets and also true in the case of datasets with more 
#'   features than samples. 
#'
#' @param data can be any dataframe or matrix with samples by column
#' @param norm.method can be prop, clr, RLE, TMM, TMMwsp
#' @param zero.remove is a value. Filter data to remove features that are 0 
#'   across at least that proportion of samples: default 0.95
#' @param zero.method can be any of NULL, prior, GBM or CZM. NULL will not 
#'   impute or change 0 values, GBM (preferred) and CZM are from the 
#'   zCompositions R package, and prior will simply add 0.5 to all counts.
#' @param log is a logical. log transform the RLE or TMM outputs, default=FALSE
#' @param group is a vector containing group information. Required for clr, RLE, 
#'   TMM, lvha, and iqlr based normalizations.
#'
#' @return Returns a list with a yes/no binary decision in 
#'   \code{is.singular} and the covariance matrix in \code{cov.matrix}
#'
#' @author Greg Gloor
#'
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom stats density
#'
#' @examples
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.singular(selex, group=group, norm.method='clr', zero.method='prior')
#' @export
aIc.singular <- function(data, norm.method='prop', zero.remove=0.95, zero.method='prior', log=FALSE, group=NULL){
  
  # remove features with 0 counts across >95% of samples 
  	data <- remove_0(data, zero.remove)
  
  # zero subustitution
  data <- zero.sub(data, zero.method)

  # aIc.get.data() is the normalization function

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
