#' \code{aIc.perturb} calculates the perturbation invariance of distance for
#'   samples with a given correction. This compares the distances of samples
#'   of the full dataset and a the perturbed dataset.
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
#' @return Returns a list with the maximum proportional perturbation in \code{ol} 
#'   (expect 0, but values up to 1% are acceptable), a yes/no binary decision in 
#'   \code{is.perturb}, the table of distances for the whole and perturbaton
#'   in \code{dist.all} and \code{dist.perturb}, the histogram of the
#'   perturbations in \code{plot}, and the plot and axis
#'   labels in \code{main} \code{xlab} and \code{ylab}. .
#'
#' @author Greg Gloor
#'
#' @export aIc.perturb
#' @importFrom zCompositions cmultRepl
#'
#' @examples
#' library(ALDEx2)
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.perturb(selex, norm.method='prop')
#' plot(x$plot, main=x$main, ylab=x$ylab, xlab=x$xlab)
aIc.perturb <- function(data, norm.method='prop', zero.method='remove', log=FALSE, group=NULL){
  
  # remove features with 0 counts across >95% of samples 
  if(zero.method == 'remove'){
  	data <- remove_0(data)
  }

  perturb <- as.vector(apply(data, 1, min) > 0)
  perturb[perturb == T] <- 5
  perturb[perturb == F] <- 1
  data.perturb <- data * perturb
  
# leave this for another day
# can implement if zero method is universally 'remove'
#  library(zCompositions)
#  data <- t(cmultRepl(t(data), method='GBM', output='p-counts'))
#  data.perturb <- t(cmultRepl(t(data.perturb), method='GBM', output='p-counts'))
  
  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  x.2 <- aIc.get.data(data.perturb, group=group, norm.method=norm.method, log=log)

  dist.all <- round(dist(t(x.1)), 4)
  dist.perturb <- round(dist(t(x.2)),4)
  
    plot.out <- hist((dist.all-dist.perturb)/dist.all *100, breaks=99, plot=F) 
    xlab='% deviance from no perturbation'
    ylab='Frequency'
    
    ol <- max(abs(plot.out$breaks))
  if(ol > 2) { 
    is.dom = 'No'
  } else { 
    is.dom = 'Yes'
  }
    main=paste('maximum absolute pertubation: ', round(ol,1),  sep="")

  return( list(ol=ol,is.perturb=is.dom, dist.all = dist.all, dist.perturb = dist.perturb, plot=plot.out,main=main, xlab=xlab, ylab=ylab))
}
