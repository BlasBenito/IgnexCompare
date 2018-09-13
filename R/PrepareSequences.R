#' Prepare input sequences.
#'
#' This function takes a pollen sequence as reference (sequence.A), and a target pollen sequence (sequence.B), checks the integrity of the reference sequence, and tries to adapt the target sequence to the characteristics of the reference sequence: matching column names and number of columns, and no missing data. Missing data can be handled in two ways: 1) deleting rows with NA or empty cases; 2) interpolation of empty data using the average or weighted average of contiguous cases.
#'
#' @param sequence.A A dataframe containing pollen counts or percentages. This sequence will be used as reference. It will be checked for empty cases.
#' @param sequence.A.name A character string with the name of the sequence.
#' @param sequence.B A dataframe containing pollen data. This sequence will be compared with and adapted to the structure of sequence.A
#' @param sequence.B.name A character string with the name of the sequence.
#' @param if.empty.cases A character argument with two possible values "omit", or "zero". The default value is "omit", but it removes every row with at least one empty record. The option "zero" replaces NA data with zeros.
#' @param transformation "none", "percentage", "proportion", "hellinger"
#' @param silent Boolean, set to TRUE to hide all messages, and set to FALSE otherwise.
#' @return A list with four slots:
#' \itemize{
#' \item \emph{taxa} Common column names of the sequences, listing the taxa included in them.
#' \item \emph{metadata} A dataframe with details about the processing of the sequences.
#' \item \emph{sequence.A} Dataframe, a processed pollen sequence.
#' \item \emph{sequence.B} Dataframe, a processed pollen sequence.
#' }
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' data(InputDataExample)
#' @export
PrepareSequences=function(sequence.A=NULL, sequence.A.name=NULL, sequence.B=NULL, sequence.B.name=NULL, if.empty.cases=NULL, transformation=NULL, silent=NULL){


  #CHECKING if.empty.cases
  ##############################################################
  if (is.null(if.empty.cases)){
    message("Argument 'if.empty.cases' is empty, setting its value to 'none'.")
    if.empty.cases="omit"
    }

  if (if.empty.cases %not-in% c("omit", "zero")){
    stop("Wrong value in the argument 'if.empty.cases'. It can only have the values 'omit' or 'zero'.")
  }


  #CHECKING transformation
  ##############################################################
  if (is.null(transformation)){
    message("Argument 'transformation' is empty, setting its value to 'none'.")
    transformation="none"
  }

  if (transformation %not-in% c("none", "percentage", "proportion", "hellinger")){
    stop("The 'transformation' argument only accepts the values: 'none', 'percentage', or 'hellinger'.")
  }

  #CHECKING INPUT DATA
  ##############################################################
  if (is.null(sequence.A)){
    stop("ERROR: reference sequence was not provided.")
  }

  if (is.null(sequence.B)){
    stop("ERROR: target sequence was not provided.")
  }

  #data not provided
  if (is.null(sequence.A.name)){
    sequence.A.name=""
    message("WARNING: reference sequence name was not provided.")
  }

  #data not provided
  if (is.null(sequence.B.name)){
    sequence.B.name=""
    message("WARNING: target sequence name was not provided.")
  }


  #SILENT?
  ##############################################################
  if (is.null(silent)){silent=FALSE}


  #TESTING DATASETS AND SUBSETTING COLUMNS
  ##############################################################
  ##############################################################
  if (silent==FALSE){cat("Checking input datasets...", sep="\n")}


  #CHECKING DATA FORMAT
  #if sequence.A is not dataframe
  if(is.data.frame(sequence.A)==FALSE){
    stop("The reference sequence is not a data frame.")
  }#end of REFERENCE SEQUENCE IS NOT DATAFRAME


  #if sequence.B is not dataframe
  if(is.data.frame(sequence.B)==FALSE){
    stop("The target sequence is not a data frame.")
  }#end of IF TARGET SEQUENCE IS NOT DATAFRAME


  #OVERLAP OF COMMON COLUMN NAMES
  common.column.names=intersect(colnames(sequence.A), colnames(sequence.B))


  #ORIGINAL DIMENSIONS OF THE DATAFRAMES
  sequence.A[is.na(sequence.A=="")]=NA
  sequence.B[is.na(sequence.B=="")]=NA
  original.na.sequence.A=sum(is.na(sequence.A))
  original.na.sequence.B=sum(is.na(sequence.B))
  original.nrow.sequence.A=nrow(sequence.A)
  original.nrow.sequence.B=nrow(sequence.B)
  original.ncol.sequence.A=ncol(sequence.A)
  original.ncol.sequence.B=ncol(sequence.B)


  #WHAT COLUMNS WERE REMOVED FROM THE TARGET DATASET?
  removed.column.names.sequence.A=setdiff(colnames(sequence.A), common.column.names)
  removed.column.names.sequence.B=setdiff(colnames(sequence.B), common.column.names)


  if (length(removed.column.names.sequence.A)==0){
    removed.column.names.sequence.A=""
  }

  if (length(removed.column.names.sequence.B)==0){
    removed.column.names.sequence.B=""
  }


  #SUBSET BY COMMON COLUMN NAMES
  sequence.A=sequence.A[, common.column.names]
  sequence.B=sequence.B[, common.column.names]


  #messages
  if (silent == FALSE){
    if (length(removed.column.names.sequence.A)==1){
      message(paste("WARNING: the column", removed.column.names.sequence.A, "was removed from the sequence A.", sep=" "))
    }

    if (length(removed.column.names.sequence.A)>1){
      message(paste("WARNING: the columns", removed.column.names.sequence.A, "were removed from the sequence A.", sep=" "))
    }

    if (length(removed.column.names.sequence.B)==1){
      message(paste("WARNING: the column", removed.column.names.sequence.B, "was removed from the sequence B.", sep=" "))
    }

    if (length(removed.column.names.sequence.B)>1){
      message(paste("WARNING: the columns", removed.column.names.sequence.B, "were removed from the sequence B.", sep=" "))
    }
  }


  #HANDLING NA DATA
  ##############################################################
  ##############################################################
  sequence.A=HandleNACases(sequence=sequence.A, sequence.name=sequence.A.name, if.empty.cases=if.empty.cases, silent=silent)
  sequence.B=HandleNACases(sequence=sequence.B, sequence.name=sequence.B.name, if.empty.cases=if.empty.cases, silent=silent)

  #counting NAs again for the metadata table
  final.na.sequence.A=sum(is.na(sequence.A))
  final.na.sequence.B=sum(is.na(sequence.B))

  #APPLYING TRANSFORMATIONS "none", "percentage", "proportion", "hellinger"
  ##############################################################
  ##############################################################

  #COMPUTING PROPORTION
  #############################
  if (transformation=="proportion"){

    sequence.A=sweep(sequence.A, 1, rowSums(sequence.A), FUN="/")
    sequence.B=sweep(sequence.B, 1, rowSums(sequence.B), FUN="/")

  }

  #COMPUTING PERCENTAGE
  ############################
  if (transformation=="percentage"){

    sequence.A=sweep(sequence.A, 1, rowSums(sequence.A), FUN="/")*100
    sequence.B=sweep(sequence.B, 1, rowSums(sequence.B), FUN="/")*100

  }

  #COMPUTING HELLINGER TRANSFORMATION
  #############################
  if (transformation=="hellinger"){

    sequence.A=sqrt(sweep(sequence.A, 1, rowSums(sequence.A), FUN="/"))
    sequence.B=sqrt(sweep(sequence.B, 1, rowSums(sequence.B), FUN="/"))

  }


  #WRAPPING UP RESULT
  ###################
  #CREATING ROWNAMES FOR DATAFRAMES (to be used for the matrixes also)
  rownames.sequence.A=1:nrow(sequence.A)
  rownames.sequence.B=1:nrow(sequence.B)

  #COERCING INTO MATRIX
  sequence.A=as.matrix(sequence.A)
  sequence.B=as.matrix(sequence.B)

  #NAMES TO COLS AND ROWS
  colnames(sequence.A)=common.column.names
  colnames(sequence.B)=common.column.names
  rownames(sequence.A)=rownames.sequence.A
  rownames(sequence.B)=rownames.sequence.B


  #PREPARING RESULTS
  ##################
  #filling metadata
  metadata=data.frame(detail=character(10),sequence.A=character(10), sequence.B=character(10), stringsAsFactors = FALSE)
  metadata[1,1:3]=c("name", sequence.A.name, sequence.B.name)
  metadata[2,1:3]=c("initial.rows", original.nrow.sequence.A, original.nrow.sequence.B)
  metadata[3,1:3]=c("final.rows", nrow(sequence.A), nrow(sequence.B))
  metadata[4,1:3]=c("initial.columns", original.ncol.sequence.A, original.ncol.sequence.B)
  metadata[5,1:3]=c("final.columns", ncol(sequence.A), ncol(sequence.B))
  metadata[6,1:3]=c("excluded.columns", removed.column.names.sequence.A, removed.column.names.sequence.B)
  metadata[7,1:3]=c("initial.empty.cases", original.na.sequence.A, original.na.sequence.B)
  metadata[8,1:3]=c("if.empty.cases", if.empty.cases, if.empty.cases)
  metadata[9,1:3]=c("final.empty.cases", final.na.sequence.A, final.na.sequence.B)
  metadata[10,1:3]=c("transformation", transformation, transformation)

  #list
  result=list()
  result[[1]]=common.column.names
  result[[2]]=metadata
  result[[3]]=sequence.A
  result[[4]]=sequence.B
  result[[5]]=NA
  result[[6]]=NA
  result[[7]]=NA
  result[[8]]=NA
  result[[9]]=NA

  names(result)=c("taxa", "metadata", "sequence.A", "sequence.B", "distance.matrix", "sum.distances.sequence.A", "sum.distances.sequence.B", "psi", "p.value")

  return(result)

}

