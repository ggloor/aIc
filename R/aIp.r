#' \code{aIc.perturb} calculates the perturbation invariance of distance for
#'   samples with a given correction. This compares the distances of samples
#'   of the full dataset and a the perturbed dataset.
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
#' @param distance can be euclidian, bray, or jaccard. euclidian on log-ratio
#'   transformed data is the same as the Aitchison distance. default=euclidian
#' @param log is a logical. log transform the RLE or TMM outputs, default=FALSE
#' @param group is a vector containing group information. Required for clr, RLE, 
#'   TMM, lvha, and iqlr based normalizations.
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
#' @examples
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' x <- aIc.perturb(selex, group=group, norm.method='clr', distance='euclidian', zero.method='prior')
#' plot(x$plot, main=x$main, ylab=x$ylab, xlab=x$xlab)
#' @export
aIc.perturb <- function(data, norm.method='prop', zero.remove=0.95, zero.method='prior', distance='euclidian', log=FALSE, group=NULL){
  
  # remove features with 0 counts across >95% of samples 
  data <- remove_0(data, zero.remove)
  	
  # zero subustitution
  data <- zero.sub(data, zero.method)

  # perturb high count values
  perturb <- as.vector(apply(data, 1, min) > 10)
  perturb[perturb == T] <- 5
  perturb[perturb == F] <- 1
  data.perturb <- data * perturb
    
  x.1 <- aIc.get.data(data, group=group, norm.method=norm.method, log=log)
  x.2 <- aIc.get.data(data.perturb, group=group, norm.method=norm.method, log=log)

  dist.all <- aIc.get.dist(x.1, distance)
  dist.perturb <- aIc.get.dist(x.2, distance)
  
    plot.out <- hist((dist.all-dist.perturb)/dist.all, breaks=99, plot=F) 
    density.out <- density((dist.all-dist.perturb)/dist.all) 
    xlab='(ref - perturb)/ref'
          
    ylab='Frequency'
    
    ol <- max(abs(plot.out$breaks))
  if(ol > 0.01) { 
    is.dom = 'No'
  } else { 
    is.dom = 'Yes'
  }
    main=paste('maximum fold change: ', round(ol,1),  sep="")

  return( list(ol=ol,is.perturb=is.dom, dist.all = dist.all, dist.perturb = dist.perturb, plot=plot.out, density=density.out, main=main, xlab=xlab, ylab=ylab))
}
