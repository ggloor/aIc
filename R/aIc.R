#' Calculate the subcompositional coherence of samples in
#' a dataset for a given correction. 
#'
#' `aIc.coherent` compares the correlation coefficients
#' of features in common of the full dataset and a subset of the dataset.
#' This is expected to be false for all compositional datasets and transforms. 
#'
#' @param data can be any dataframe or matrix with samples by column
#' @param norm.method can be prop, clr, RLE, TMM, TMMwsp, lvha, iqlr
#' @param zero.remove is a value. Filter data to remove features that are 0 
#' across at least that proportion of samples: default 0.95
#' @param zero.method can be any of NULL, prior, GBM or CZM. NULL will not 
#' impute or change 0 values, GBM and CZM are from the 
#' zCompositions R package, and prior will simply add 0.5 to all counts.
#' @param log is a logical. log transform the prop, RLE or TMM outputs, default=FALSE
#' @param group is a vector containing group information. Required for clr, RLE, 
#' @param cor.test is either the pearson or spearman method (default)
#'
#' @return Returns a list with the correlation in \code{cor}, a yes/no binary 
#' decision in \code{is.coherent},  the x and y values for a scatterplot
#' of the correlations in the full and subcompositions, and the plot and axis
#' labels in \code{main} \code{xlab} and \code{ylab}. 
#'
#' @author Greg Gloor
#'
#' @importFrom grDevices rgb
#' @importFrom graphics abline hist
#' @importFrom stats cor cov dist runif
#' @importFrom zCompositions cmultRepl lrSVD
#' @importFrom vegan vegdist
#' @examples
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.coherent(selex, group=group, norm.method='clr', zero.method='prior')
#' plot(x$plot[,1], x$plot[,2], main=x$main, ylab=x$ylab, xlab=x$xlab)
#'
#' @export
aIc.coherent <- function(data, norm.method="prop", zero.remove=0.95, zero.method='prior', log=FALSE, group=NULL, cor.test='spearman'){

  # remove features with 0 counts across >95% of samples 
  data <- remove_0(data, zero.remove)

  # zero subustitution
  data <- zero.sub(data, zero.method)

  # aIc.get.data() is the normalization function
 
  size.sub <- floor(nrow(data)/2)
  data.sub <- data[1:size.sub,]

  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  x.2 <- aIc.get.data(data.sub, group=group, norm.method=norm.method, log=log)
  
  c.x1 <- stats::cor(t(x.1), method=cor.test)
  c.x2 <- stats::cor(t(x.2), method=cor.test)
  
  v.x1 <- c(c.x1[1:size.sub,1:size.sub])
  v.x2 <- c(c.x2)
  
  # convert the correlations to a vector
  cohere <- cor(v.x1, v.x2)
  is.comp <- vector()
  if(cohere < 1){ 
    is.comp = 'No'
  } else {
    is.comp = 'Yes'
  }
  
  r.num = floor(runif(10000, min=1, max=length(v.x1)))
  
  if( length(v.x1) < length(r.num) ) {
    plot.out <- cbind(v.x1, v.x2)
  } else {
    plot.out <- cbind(v.x1[r.num], v.x2[r.num])
  }	
  
  colnames(plot.out) <- c('sub', 'full')
  main=paste('correlation between full and sub is: ', round(cohere,4))
  xlab='sub composition'
  ylab='full composition' 
  
  return(list(cor = cohere, is.coherent = is.comp, plot=plot.out, main=main, xlab=xlab, ylab=ylab))
}

