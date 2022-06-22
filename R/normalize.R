#' @importFrom edgeR DGEList calcNormFactors
# convenience functions

aIc.get.data <- function(data, norm.method=norm.method, group=NULL, log=FALSE){

  valid.norm = c('none', 'prop', 'clr', 'TMM', 'TMMwsp', 'RLE')
  if(!(norm.method %in% valid.norm)) stop('valid normalizations are none, clr, prop, TMM, TMMwsp, RLE')
 
  if(norm.method == 'none'){
    if(log==F) { x <- as.matrix(data) }
    if(log==T) { stop ('cannot take logarithm of 0 values') }
  }
 if(norm.method == 'clr'){
    x <- aIc.get.norm(data, norm.method='clr')
  } else if(norm.method == 'prop'){
    x <- aIc.get.norm(data, norm.method='prop', log=log)
  } else if(norm.method %in% c('TMM', 'RLE', 'upperquartile', 'TMMwsp')){
    x <- aIc.get.norm(data, norm.method, group = group, log=log)
  }
  if(any(is.na(x))) stop('too many 0 values. Consider adding a small pseudocount prior to running')
  if(any(is.infinite(x))) stop('too many 0 values. Consider adding a small pseudocount prior to running')
  return(x)
}

aIc.get.norm <- function(data, norm.method, group, log=FALSE){	
  
	if ( norm.method == 'clr' ){
		return(apply(data+0.5, 2, function(x) log(x) - mean(log(x))))
	}
	else if( norm.method == 'prop'){
		if(log==F){
		  return(apply(data +0.5, 2, function(x) x/sum(x)))	
		}else{
		  return( apply(data +0.5, 2, function(x) log( x/sum(x) ) )	)
		}
	}
	else {
		if(log==F){
	  	y <- edgeR::DGEList(counts=data, group=group)
	  	y <- edgeR::calcNormFactors(y, method=norm.method)
			return(t(apply(y[[1]], 1, function(x) x * y[[2]]$norm.factors)))
		} else {
	  	y <- edgeR::DGEList(counts=data+0.5, group=group)
	  	y <- edgeR::calcNormFactors(y, method=norm.method)
			return(log(t(apply(y[[1]], 1, function(x) x * y[[2]]$norm.factors))))
		}
	}
}
