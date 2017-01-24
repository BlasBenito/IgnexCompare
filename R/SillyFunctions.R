#' Small silly functions
#'
#' Two small helper functions: \strong{%not-in%}, as to check "is 1 in [2, 3, 4]?", and \strong{FirstUp}, used to turn the first letter of a genus name into uppercase.
#'
#' @param rows Number of rows of the results table.
#' @return A dataframe with n rows, and a predefined set of columns
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #example for %not-in%
#' if (1 %not-in% c(2, 3, 4)) print("1 is not in (2, 3, 4)")
#'
#' #example for FirstUp
#' FirstUp("lynx)
#'
#' @export
'%not-in%' = function(x,y)!('%in%'(x,y))

#' @export
FirstUp = function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
