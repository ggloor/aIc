#' \code{aIc.scale} calculates the scaling invariance of a sample in
#'   a dataset for a given correction. This compares the distances of samples
#'   of the full dataset and a scaled version of the dataset.
#'   This is expected to be true if the transform is behaving rationally in
#'   compositional datasets. 
#'
#' @param data can be any dataframe or matrix with samples by column
#' @param norm.method can be prop, clr, RLE, TMM, TMMwsp
#' @param zero.method is a logical. Filter data to remove features that are 0
#'   across all samples. Default is TRUE
#' @param log is a logical. log transform the RLE or TMM outputs, default=FALSE
#' @param group is a vector containing group information. Required fro RLE and 
#'   TMM based normalizations.
#'
#' @return Returns a list with the overlap between distances in the full and 
#'   scaled composition in \code{ol} (expect 0), a yes/no binary decision in 
#'   \code{is.scale} and the table of distances for the whole and scaled composition
#'   in \code{dist.all} and \code{dist.scale}, a plot showing a histogram of the resulting
#'   overlap in distances in \code{plot}, and the plot and axis
#'   labels in \code{main} \code{xlab} and \code{ylab}
#'
#' @author Greg Gloor
#'
#' @export aIc.scale
#'
#' @examples
#' library(ALDEx2)
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.scale(selex, norm.method='prop')
#' plot(x$plot, main=x$main, ylab=x$ylab, xlab=x$xlab)
aIc.scale <- function(data, norm.method='prop', zero.method='remove', log=FALSE, group=NULL){
  
  # remove features with 0 counts across >95% of samples 
  if(zero.method == 'remove'){
  	data <- remove_0(data)
  }
  
  if( min(data) == 0) {
    data <- t(zCompositions::cmultRepl(t(data), method='GBM', output='p-counts', z.warning=.99, suppress.print=TRUE))
  }  
  
  # scale by 5-fold change
  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  x.2 <- aIc.get.data(data * 5, group=group, norm.method=norm.method, log=log)

  dist.all <- dist(t(x.1))
  dist.scale <- dist(t(x.2))
  
  plot.out <- hist(dist.all-dist.scale, breaks=99, plot=F) #, 
  xlab='difference between normal and scaled composition distance'
  ylab='Frequency'

    ol <- max(abs(plot.out$breaks))
  if(ol > 2) { 
    is.dom = 'No'
  } else { 
    is.dom = 'Yes'
  }
    main=paste('maximum absolute pertubation: ', round(ol,1),  sep="")

  return( list(ol=ol,is.scale=is.dom, dist.all = dist.all, dist.scale = dist.scale, plot=plot.out, main=main, xlab=xlab, ylab=ylab))
}
