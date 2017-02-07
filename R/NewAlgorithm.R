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

distance.matrix=sequences$distance.matrix

#number of cells from corner to corner
n.cells=sum(dim(distance.matrix))-1

a=vector() #rows
b=vector() #columns
cost=vector()

#convert to dataframe
for (column in 1:ncol(distance.matrix)){
  for (row in 1:nrow(distance.matrix)){
    a=c(a, row)
    b=c(b, column)
    cost=c(cost, distance.matrix[row, column])
  }
}

dm.df=data.frame(a=a, b=b, cost=cost)

#sorting dataframe
dm.df=dm.df[order(dm.df$a, dm.df$b, dm.df$cost), ]

#removing cases that do not follow the moving rules
