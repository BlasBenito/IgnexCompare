#' Scales pollen data
#'
#' Scales pollen data between a minimum and a maximum value.
#'
#' @param sequence A dataframe containing pollen data. This sequence will be compared with and adapted to the structure of sequence.A
#' @param columns Optional, character vector of column names. If not provided, the function will try with all columns.
#' @param new.min Minimum new value, 0 by default.
#' @param new.max Maximum new value, 100 by default.
#' @return The function returns the reformatted target sequence as a data frame.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' data(InputDataExample)
#' @export
#'
RescalePollenSequence=function(sequence, columns=NULL, new.min=NULL, new.max=NULL){

  #default value for if.empty.cases
  if (is.null(sequence)){
    stop("No sequence was provided.")
  }

  if (is.null(columns)){
    columns=colnames(sequence)
  }

  if (is.null(new.min)){
    new.min=0
  }

  if (is.null(new.max)){
    new.max=100
  }

    #iterating through matching columns
  for (i in 1:length(columns)){

    column=columns[i]

    #data extremes
    old.min=min(sequence[, column])
    old.max=max(sequence[, column])

    #avoid dividing by 0
    if (old.max==0){old.max=0.0001}

    #scaling
    sequence[, column]=((sequence[, column] - old.min) / (old.max - old.min)) * (new.max - new.min) + new.min

  }

return(sequence)

}
