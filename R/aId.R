#' \code{aIc.dominant} calculates the subcompositional dominance of a sample in
#'   a dataset for a given correction. This compares the distances of samples
#'   of the full dataset and a subset of the dataset.
#'   This is expected to be true if the transform is behaving rationally in
#'   compositional datasets. 
#'
#' @param data can be any dataframe or matrix with samples by column
#' @param norm.method can be prop, clr, RLE, TMM, TMMwsp
#' @param zero.remove is a logical. Filter data to remove features that are 0 
#'   across a proportion of samples over 0.95. Default=TRUE
#' @param zero.method can be any of NULL, prior, GBM or CZM. NULL will not 
#'   impute or change 0 values, GBM (preferred) and CZM are from the 
#'   zCompositions R package, and prior will simply add 0.5 to all counts.
#' @param log is a logical. log transform the RLE or TMM outputs, default=FALSE
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
#' library(ALDEx2)
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.dominant(selex, group=group, norm.method='clr', zero.method='prior')
#' plot(x$plot, main=x$main, ylab=x$ylab, xlab=x$xlab)
#' @export
aIc.dominant <- function(data, norm.method='prop', zero.remove=TRUE, zero.method=NULL, 
  log=FALSE, group=NULL){
  
  # remove features with 0 counts across >95% of samples 
  if(zero.remove == TRUE){
  	data <- remove_0(data)
  }
  # zero subustitution
  data <- zero.sub(data, zero.method)

  # aIc.get.data() is the normalization function

  size.sub <- floor(nrow(data)/2)
  data.sub <- data[1:size.sub,]
  
  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  x.2 <- aIc.get.data(data.sub, group=group, norm.method=norm.method, log=log)

  dist.all <- dist(t(x.1))
  dist.sub <- dist(t(x.2))
  
  ol <- 100 * min(c(sum(dist.all-dist.sub < 0)/length(dist.sub),sum(dist.all-dist.sub > 0)/length(dist.sub) ))
  if(ol > 0) { 
    is.dom = 'No'
  } else { 
    is.dom = 'Yes'
  }
  plot.out <- hist((dist.all-dist.sub)/dist.all * 100, breaks=99, plot=F) #, 
  main=paste('% of dominant distances ', 100 - round(ol, 3), sep="")
  xlab='difference between full and sub composition distance'
  ylab='Frequency'

  return( list(ol=ol,is.dominant=is.dom, dist.all = dist.all, dist.sub = dist.sub, plot=plot.out, main=main, xlab=xlab, ylab=ylab))
}
