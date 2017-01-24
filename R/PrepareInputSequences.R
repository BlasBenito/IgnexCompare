#' Prepare input sequences.
#'
#' This function takes a pollen sequence as reference (reference.sequence), and a target pollen sequence (target.sequence), checks the integrity of the reference sequence, and tries to adapt the target sequence to the characteristics of the reference sequence: matching column names and number of columns, and no missing data. Missing data can be handled in two ways: 1) deleting rows with NA or empty cases; 2) interpolation of empty data using the average or weighted average of contiguous cases.
#'
#' @param reference.sequence A dataframe containing pollen counts or percentages. This sequence will be used as reference. It will be checked for empty cases.
#' @param target.sequence A dataframe containing pollen data. This sequence will be compared with and adapted to the structure of reference.sequence
#' @param if.empty.cases A character argument with two possible values "omit" or "interpolate".
#' @param taxon.agnostic A boolean argument. TRUE if the user doesn't care about the column names (pollen types) being equal between sequences, and FALSE otherwise. The default is FALSE.
#' @return The function returns the reformatted target sequence as a data frame.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' data(InputDataExample)

#' @export
PrepareInputSequences=function(reference.sequence=NULL, target.sequence=NULL, if.empty.cases=NULL, taxon.agnostic=NULL){

  stop.status=FALSE

  #defining the default value for fuzzy match
  if (is.null(taxon.agnostic)){taxon.agnostic=FALSE}
  if (is.null(if.empty.cases)){if.empty.cases="omit"}

  #data not provided
  if (is.null(reference.sequence)){
    stop("FATAL ERROR: reference sequence not provided.")
  }

  if (is.null(target.sequence)){
    stop("FATAL ERROR: target sequence not provided.")
  }



  #TESTING DATASETS
  cat("Checking input datasets...")


  #REFERENCE SEQUENCE IS NOT DATAFRAME
  if(is.data.frame(reference.sequence)==FALSE){
    message("ERROR: The reference sequence is not a data frame.")
    stop.status=TRUE
  }#end of REFERENCE SEQUENCE IS NOT DATAFRAME


  #IF TARGET SEQUENCE IS NOT DATAFRAME
  if(is.data.frame(target.sequence)==FALSE){
    message("ERROR: The target sequence is not a data frame.")
    stop.status=TRUE
  }#end of IF TARGET SEQUENCE IS NOT DATAFRAME


  #IF REFERENCE SEQUENCE IS DATAFRAME
  if(is.data.frame(reference.sequence)==TRUE){
    cat("Check correct: the reference sequence is a data frame.", sep="\n")
  }#end of IF REFERENCE SEQUENCE IS DATAFRAME


    #IF TARGET SEQUENCE IS DATAFRAME
    if(is.data.frame(target.sequence)==TRUE){
      cat("Check correct: the target sequence is a data frame.", sep="\n")
    }#end of IF TARGET SEQUENCE IS DATAFRAME


  #NOT TAXON AGNOSTIC
  if (taxon.agnostic==FALSE){

    #CHECKING COLUMN NAMES
    if (ncol(reference.sequence)==ncol(target.sequence)){
      cat("Check correct: both sequences have the same number of columns.", sep="\n")

        #IF COLNAMES ARE THE SAME
        if (colnames(reference.sequence)==colnames(target.sequence)){
          cat("Check correct: both sequences have the same column names.", sep="\n")
        }#end of IF COLNAMES ARE THE SAME


        #IF COLNAMES ARE DIFFERENT
        if (colnames(reference.sequence)!=colnames(target.sequence)){

          #overlap of column names
          common.column.names=intersect(colnames(reference.sequence), colnames(target.sequence))

          #SUBSET TARGET DATASET
          if(common.column.names==colnames(reference.sequence) & length(common.column.names)< length(colnames(target.sequence)))){

            #WHAT COLUMNS WERE REMOVED FROM THE TARGET DATASET?
            removed.colnames=colnames(target.sequence) %not-in% common.column.names
            removed.colnames=colnames(target.sequence)[removed.colnames]

            #messages
            if (length(removed.colnames)==1){
              message(paste("WARNING: the column", removed.colnames, "was removed from the target dataset", sep=" "))
            }

            if (length(removed.colnames)>1){
              message(paste("WARNING: the following columns were removed from the target dataset:", removed.colnames, sep=" "))
            }#end of messages

            #subset on target sequence
            target.sequence=target.sequence[, common.column.names]

          }#end of SUBSET TARGET DATASET


        }#end of IF COLNAMES ARE DIFFERENT

    }#end of CHECKING COLUMN NAMES

  }#end of NOT TAXON AGNOSTIC




  if (stop.status==TRUE){
    stop("ERROR: Data errors could not be handled. Please, reformat your data and try again. You have two lives left.")
  }


  if (stop.status==FALSE){
    cat("Returning processed target sequence.", sep="\n")
    return(target.sequence)
  }


}
