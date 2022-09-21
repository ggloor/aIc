#' \code{aIc.dominant} calculates the subcompositional dominance of a sample in
#'   a dataset for a given correction. This compares the distances of samples
#'   of the full dataset and a subset of the dataset.
#'   This is expected to be true if the transform is behaving rationally in
#'   compositional datasets. 
#'
#' @param data can be any dataframe or matrix with samples by column
#' @param norm.method can be prop, clr, RLE, TMM, TMMwsp
#' @param zero.remove is a value. Filter data to remove features that are 0 
#'   across at least that proportion of samples: default 0.95
#' @param zero.method can be any of NULL, prior, GBM or CZM. NULL will not 
#'   impute or change 0 values, GBM (preferred) and CZM are from the 
#'   zCompositions R package, and prior will simply add 0.5 to all counts.
#' @param log is a logical. log transform the RLE or TMM outputs, default=FALSE
#' @param distance can be euclidian, bray, or jaccard. euclidian on log-ratio
#'   transformed data is the same as the Aitchison distance. default=euclidian
#' @param group is a vector containing group information. Required for clr, RLE, 
#'   TMM, lvha, and iqlr based normalizations.
#'
#' @return Returns a list with the overlap between distances in the full and 
#'   subcompositon in \code{ol} (expect 0), a yes/no binary decision in 
#'   \code{is.dominant} and the table of distances for the whole and subcomposition
#'   in \code{dist.all} and \code{dist.sub}, a plot showing a histogram of the resulting
#'   overlap in distances in \code{plot}, and the plot and axis
#'   labels in \code{main} \code{xlab} and \code{ylab}
#'
#' @author Greg Gloor
#'
#' @examples
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.dominant(selex, group=group, norm.method='clr', distance='euclidian', zero.method='prior')
#' plot(x$plot, main=x$main, ylab=x$ylab, xlab=x$xlab)
#' @export
aIc.dominant <- function(data, norm.method='prop', zero.remove=0.95, zero.method='prior', 
  log=FALSE, distance='euclidian', group=NULL){
  
  # remove features with 0 counts across >95% of samples 
  data <- remove_0(data, zero.remove)
  
  # zero substitution
  data <- zero.sub(data, zero.method)

  # aIc.get.data() is the normalization function

  size.sub <- floor(nrow(data)/2)
  data.sub <- data[1:size.sub,]
  
  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  x.2 <- aIc.get.data(data.sub, group=group, norm.method=norm.method, log=log)

  dist.all <- aIc.get.dist(x.1, distance)
  dist.sub <- aIc.get.dist(x.2, distance)
    
#  ol <- min(c(sum(dist.all-dist.sub < 0)/length(dist.sub),sum(dist.all-dist.sub > 0)/length(dist.sub) ))
  ol <- 1 - (sum((dist.all-dist.sub)/dist.all <0) / length(dist.all))
  
  if(ol < 1) { 
    is.dom = 'No'
    main=paste('Proportion of dominant distances ', round(ol, 3), sep="")
  } else { 
    is.dom = 'Yes'
    main=paste('Proportion of dominant distances ', round(ol, 3), sep="")
  }
  plot.out <- hist((dist.all-dist.sub)/dist.all, breaks=99, plot=F) #, 
  density.out <- density((dist.all-dist.sub)/dist.all) #, 
  xlab='Relative distance between full and sub composition '
  ylab='Frequency'

  return( list(ol=ol,is.dominant=is.dom, dist.all = dist.all, dist.sub = dist.sub, plot=plot.out, density=density.out, main=main, xlab=xlab, ylab=ylab))
}
