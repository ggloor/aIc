#' \code{aIc.plot} plots the result of the distance tests.
#'
#' @param test.out is the output from either aIc.dominant, aIc.scale, aIc.perturb
#'
#' @return returns a plot of the density of the distance test results.
#    This includes summary information and a graphical interpretation of the 
#'   test result.
#'
#' @author Greg Gloor
#' 
#' @importFrom graphics text abline rect layout par
#'
#' @examples
#' data(selex)
#' group = c(rep('N', 7), rep('S', 7))
#' test.out <- aIc.dominant(selex, norm.method='prop', group=group)
#' aIc.plot(test.out)
#' @export
aIc.plot <- function(test.out){
  x <- test.out
  
  if(names(x)[2] == "is.dominant"){
     min.dens <- min(x$density$x)
     if(min.dens > 0){min.dens = 0}
     max.dens <- max(x$density$x)
     if(max.dens < 0){max.dens = 0}
     
     plot(x$density, xlim=c(min.dens, max.dens), main=x$main, xlab=x$xlab)
     abline(v=0, lwd=3, lty=2, col="red")
       text(max(x$density$x), max(x$density$y) *.8, adj=1, 
         labels="dominant", col='red')
       if(x$ol != 0){
         text(min(x$density$x) *.5 , max(x$density$y) *.8, adj=1, 
           labels="not dominant", col='red')
        }   
  } else if(names(x)[2] == "is.perturb" || names(x)[2] == "is.scale"){
     min.dens <- min(x$density$x)
     max.dens <- max(x$density$x)
     
     if( min.dens < 0 && max.dens < 0){max.dens = abs(min.dens) *.25}
     if( min.dens > 0 && max.dens > 0){min.dens = 0 - abs(max.dens) *.25}
     
     
     plot(x$density, xlim=c(min.dens, max.dens),main=x$main, xlab=x$xlab, 
       ylim=c(0, max(x$density$y) * 1.1))
     rect(-0.005, 0, 0.005, max(x$density$y), col=rgb(1,0,0,0.2), border=NA)
     if(x[[2]] == "No") {
       text(0, max(x$density$y) *1.05, labels='not invariant')
     } else {
       text(0, max(x$density$y) *1.05, labels='invariant')
     }
  } else if(names(x)[2] == "is.coherent"){
  
     # reset parameters on exit from this 
     opar <- par(mar=c(5.1, 4.1, 4.1, 2.1), mfrow=c(1,1))
     on.exit(par(opar))

     layout(matrix(c(2,0,1,3), nrow=2, ncol=2, byrow=T), widths=c(4,1), 
       heights=c(1,4), respect=T)

	# mar(B,L,T,R)
	# top and right of main plot to 0 margin

	par(mar=c(5.1,4.1,0,0))
	plot(x$plot, pch=19, col=rgb(0,0,0,0.2), ces=0.5, xlim=c(-1,1), 
	  ylim=c(-1,1))
    abline(0,1, lty=2, lwd=3, col='red')
     if(x[[2]] == "No") {
       text(-0.5, 0.8 *1.05, labels='not coherent')
     } else {
       text(-0.5, 0.8 *1.05, labels='coherent')
     }

	par(mar=c(0,4.1,0,0))
	d.jnk <- density(x$plot[,1])
	plot(d.jnk$x, d.jnk$y, type='l', col='grey40', xlab=NA, ylab=NA, bty='n', xaxt="n", yaxt='n', lwd=2)

	par(mar=c(5.1,0,0,0))
	d.jnk <- density(x$plot[,1])
	plot(d.jnk$y, d.jnk$x, type='l', col='grey40', xlab=NA, ylab=NA, bty='n', xaxt='n', yaxt='n', lwd=2)
    
  } else {
  	stop("this will only plot for the dominance, scale, coherence or perturbation tests")
  }
}