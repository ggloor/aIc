#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom stats as.dist
#' @importFrom ALDEx2 aldex.set.mode
# convenience functions

aIc.get.data <- function(data, norm.method=norm.method, group=NULL, log=FALSE){

  if(norm.method=='clr'){
    denom='all'
  } else if(norm.method=='iqlr') { 
    denom=norm.method
    norm.method='clr'
  } else if(norm.method=='lvha') { 
    denom=norm.method
    norm.method='clr'
  }

  valid.norm = c('none', 'prop', 'clr', 'TMM', 'TMMwsp', 'RLE', 'iqlr', 'lvha')
  if(!(norm.method %in% valid.norm)) stop('valid normalizations are none, clr, iqlr, lvha, prop, TMM, TMMwsp, RLE')
  
  if(norm.method == 'none'){
    if(log==F) { x <- as.matrix(data) }
    if(log==T) { stop ('cannot take logarithm of 0 values') }
  }
  
 if(norm.method == 'clr'){
 	  if(min(data == 0)){ stop('please choose a 0 replacement method') }
 	  if(is.null(group)){ stop('the group vector must be defined') }
    den.vec <- ALDEx2::aldex.set.mode(data, group, denom)
    x <- aIc.get.norm(data, norm.method='clr', group=group, den.vec=den.vec)
  } else if(norm.method == 'prop'){
    if(min(data)==0 & log==T) { stop ('cannot take logarithm of 0 values, choose a 0 replacement method') }
    x <- aIc.get.norm(data, norm.method='prop', log=log, den.vec=NULL)
  } else if(norm.method %in% c('TMM', 'RLE', 'upperquartile', 'TMMwsp')){
    if(min(data) == 0 & log==T) { stop ('cannot take logarithm of 0 values, choose a 0 replacement method') }
    if(is.null(group)){ stop('the group vector must be defined') }
    x <- aIc.get.norm(data, norm.method, group=group, log=log, den.vec=NULL)
  }
  if(any(is.na(x))) stop('too many 0 values. Choose a 0 replacement method')
  if(any(is.infinite(x))) stop('too many 0 values. Choose a 0 replacement method')
  return(x)
}

aIc.get.norm <- function(data, group, norm.method, log=FALSE, den.vec=NULL){	

	if ( norm.method == 'clr' ){
		if(is.null(den.vec)){stop('den.vec is null obiwan')} ##
    return(apply(data, 2, function(x) log(x) - mean(log(x[den.vec]))))
	}
	else if( norm.method == 'prop'){
		if(log==F){
		  return(apply(data, 2, function(x) x/sum(x)))	
		}else{
		  return( apply(data, 2, function(x) log( x/sum(x) ) )	)
		}
	}
	else {
	  	y <- edgeR::DGEList(counts=data, group=group)
	  	y <- edgeR::calcNormFactors(y, method=norm.method)
		if(log==F){
			return(t(apply(y[[1]], 1, function(x) x * y[[2]]$norm.factors)))
		} else {
			return(log(t(apply(y[[1]], 1, function(x) x * y[[2]]$norm.factors))))
		}
	}
}

# removes 0 only rows
# filters by at least 95% non-0 occurrence
remove_0 <- function(data, zero.remove){
    if(!is.numeric(zero.remove)) { stop( 'enter a valid proportion for zero removal' ) }
		if(zero.remove < 0 | zero.remove > 1){ stop( 'enter a valid proportion for zero removal' ) }

		data <- data[rowSums(data) > 0,]
		n0 <- apply(data, 1, function(x) length(which(x == 0)))
		max.0 = ncol(data) * zero.remove
		return(data[n0 < max.0,])      
}

# zero subsition, or not
zero.sub <- function(data, zero.method){
   if(is.null(zero.method)){
    warning('0 values are not substituted or imputed in any way')
    data <- data 
	} else if( zero.method=='GBM' & min(data) == 0) {
    data <- t(zCompositions::cmultRepl(t(data), method='GBM', output='p-counts', 
      z.warning=.99, suppress.print=TRUE))
  } else if( zero.method=='CZM' & min(data) == 0) {
    data <- t(zCompositions::cmultRepl(t(data), method='CZM', output='p-counts', 
      z.warning=0.65, suppress.print=TRUE))
  }else if( zero.method=='lrSVD' & min(data) == 0) {
  	dl <- rep(0.5, nrow(data))
    data <- t(zCompositions::lrSVD(t(data), label=0, dl=dl, 
      z.warning=0.85))
  } else if( zero.method=='prior' & min(data) == 0) {
		data <- data + 0.5
  }
  return(data)
}

# convenience distance function
# shared by aIc.dominant, aIc.perturb, aIc.scale
aIc.get.dist <- function(x, distance){
  if(distance == 'euclidian' ){
    return(dist(t(x), distance) )
  }else if(distance == 'bray' | distance == 'jaccard'){
  	return(vegan::vegdist(t(x), distance))
  } else if(distance == 'jsd'){
    # Jensen-Shannon distance
    return(dist.JSD(x))
 }
}

# actually square root of Jensen-Shannon Divergence
# from http://enterotype.embl.de/enterotypes.html

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- ncol(inMatrix)
	matrixRowSize <- nrow(inMatrix)
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) {
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	stats::as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix)
}