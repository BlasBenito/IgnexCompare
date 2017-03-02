#' running description
#'
#' description
#'
#' @usage
#'
#' @param param Number of rows of the results table.
#' @return whatever
#' \itemize{
#'  \item{"parameter 1"}{Stuff}
#'  \item{"parameter 2"}{Stuff}
#' }
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
PlotSlotting=function(slotting, main=NULL){

  #TO DO: WRITE CHECKS FOR OBJECT NAMES
  if (is.null(main)){main="Sequence slotting."}

  distance.matrix=slotting$distance.matrix
  best.distances=slotting$lowest.costs
  best.solution=slotting$best.slotting
  psi=slotting$psi
  p.value=slotting$p.value

  #PLOTTING
  par(mfrow=c(2,1), mar=c(3,4,2,2), oma=c(1,1,2,1))
  PlotDistanceMatrix(distance.matrix, main=paste("Best slotting: psi =", round(psi, 2),"; p-value = ", p.value, sep=" "), path=best.solution)
  plot(best.distances, type="l", xlab="Best solutions", ylab="Overall cost", main="Slotting improvement", lwd=2, col="red4")
  mtext(main, 3, outer=TRUE, cex=1.5)

}


