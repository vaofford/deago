#' @title Construct input matrix for DESeq2
#' @description Rows of coldata correspond to columns of countdata
#'
#' @param targets Name of object containing target data
#' @param countdata Name of object containing count data
#' @param parameters Parameters
#'
#' @importFrom stats relevel
#' @return input matrix for DESeq2 analysis
#' @export

prepareColData <- function(targets, countdata, parameters)
{
  condition <- factor( tolower(targets$condition) )

  # Relevel if ref given - will default to alphabetical (see order in importTargets)
  if ( "control" %in% names(parameters)) {
    if (!is.null(parameters$control)) {
      ref <- parameters$control
      ref <- tolower(ref)
      stopifnot(ref %in% condition)
      condition <- relevel(condition, ref=ref)
    }
  }

  coldata <- data.frame(row.names=colnames(countdata), condition)
  return(coldata)
}

#' @title Construct DESeq2 dataset
#' @description A DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values, intermediate calculations and results of an analysis of differential expression.
#'
#' @param countdata Name of object containing count data
#' @param coldata input matrix for DESeq2 analysis
#'
#' @import DESeq2
#'
#' @return Returns a DESeqDataSet object
#' @export

constructDESeqDataset <- function(countdata, coldata)
{
  # This is a seperate function as will need to change the design if there are multiple factors/interactions or batch effects included

  # Prepare DESeq2 matrix
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

  return(dds)
}

#' @title Differential expression analysis using DESeq2
#' @description Performs a default DESeq2 analysis
#'
#' @param targets Name of object containing target data
#' @param countdata Name of object containing count data
#' @param parameters Levels are re-ordered so that the level specified by ref is first
#'
#' @import DESeq2
#'
#' @return Returns a DESeqDataSet object
#' @export

runDESeqAnalysis <- function(targets, countdata, parameters)
{
  # Construct input matrix for DESeq2
  coldata <- prepareColData(targets,countdata, parameters)

  # Construct DESeq2 dataset
  dds <- constructDESeqDataset(countdata, coldata)

  # Run default analysis
  dds <- DESeq(dds)

  #return(dds)
}

