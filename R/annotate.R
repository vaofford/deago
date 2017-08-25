#' @title Annotate DESeqDataSet
#'
#' @description \code{annotateDataset} reads a tab-delimited annotation file
#'   containing gene symbols and/or GO terms into a dataframe which is used to
#'   annotate the DESeqDataSet.
#'
#' @details The annotation file is a text file containing gene information. Rows
#'   correspond to genes while columns correspond to gene symbols and GO terms
#'   that are associated with the genes. The annotation file should \emph{not}
#'   contain a header row.
#'
#'   The first column must contain the gene identifiers which correspond to the
#'   row names of the \link[DESeq2]{DESeqDataSet} object. There doesn't need to
#'   be an entry for every gene identifier in the \link[DESeq2]{DESeqDataSet}
#'   but, at least one identifier should be represented.
#'
#'   There then follows one or two further columns depending on whether the
#'   annotations consist of only gene symbols, only GO terms or both gene
#'   symbols and GO terms. Where there are multiple GO terms for a gene, the
#'   individual terms should be semi-colon (;) delimited e.g.
#'   \emph{GO:0004415;GO:0004415;GO:0004415;}.
#'
#'   \code{annotateDataset} returns a \link[DESeq2]{DESeqDataSet} with the
#'   annotations stored in the object metadata. Gene symbols will be accessible
#'   using \code{mcols(dds)$symbol} while GO terms will be accessible using
#'   \code{metadata(dds)$go} where \code{dds} represents a
#'   \link[DESeq2]{DESeqDataSet} object.
#'
#' @family annotation functions
#' @family import functions
#'
#' @param dds A \link[DESeq2]{DESeqDataSet} object.
#' @param parameters A \link{list} containing key/value pairs which define the
#'   parameters for the analysis (see \link[deago]{importConfig} for more
#'   information).
#'
#' @return A \link[DESeq2]{DESeqDataSet} object with annotations.
#'
#' @importFrom methods is
#' @importFrom DESeq2 DESeqDataSet
#'
#' @export
#'
#' @examples
#' \dontrun{
#' annotateDataset(dds, parameters)
#' }

annotateDataset <- function(dds, parameters)
{
  stopifnot(is(dds, 'DESeqDataSet'))
  stopifnot(is.list(parameters))
  if(!exists('annotation_file', where=parameters)) stop ("Could not import annotation: no annotation file provided.")

  ann <- importAnnotation(parameters$annotation_file)

  if (length(intersect(rownames(dds), ann[,1])) == 0) stop("Could not import annotation: no dataset identifiers found.")
  dds <- addAnnotationsToDataSet(dds, ann)

  return(dds)
}

#' @title Import Annotation File
#'
#' @description \code{importAnnotation} reads a tab-delimited annotation file
#'   containing gene symbols and/or GO terms into a dataframe.
#'
#' @details The annotation file is a text file containing gene information. Rows
#'   correspond to genes while columns correspond to gene symbols and GO terms
#'   that are associated with the genes. The annotation file should \emph{not}
#'   contain a header row.
#'
#'   The first column must contain the gene identifiers which correspond to the
#'   row names of the \link[DESeq2]{DESeqDataSet} object. There doesn't need to
#'   be an entry for every gene identifier in the \link[DESeq2]{DESeqDataSet}
#'   but, at least one identifier should be represented.
#'
#'   There then follows one or two further columns depending on whether the
#'   annotations consist of only gene symbols, only GO terms or both gene
#'   symbols and GO terms. Where there are multiple GO terms for a gene, the
#'   individual terms should be semi-colon (;) delimited e.g.
#'   \emph{GO:0004415;GO:0004415;GO:0004415;}.
#'
#' @family annotation functions
#' @family import functions
#'
#' @param file A \link{character} string giving the name of the annotation file.
#'
#' @return A dataframe with annotations.
#'
#' @importFrom utils read.delim
#'
#' @export
#'
#' @examples
#' \dontrun{
#' importAnnotaton(annotation_file)
#' }

importAnnotation <- function(file)
{
  stopifnot(is(file, 'character'))
  if (!file.exists(file)) stop(paste0(file," does not exist"))

  ann <- read.delim(file=file, header = FALSE, quote = "", sep ="\t", colClasses = "character")

  if (ncol(ann) < 2) stop("Could not import annotation: too few columns in annotation file.")
  if (ncol(ann) > 3) stop("Could not import annotation: too many columns in annotation file.")

  return(ann)
}

#' @title Add Annotations To DESeqDataSet
#'
#' @description \code{addAnnotationsToDataSet} adds annotations from a dataframe
#'   containing gene symbols and/or GO terms to the \link[DESeq2]{DESeqDataSet}
#'   object.
#'
#' @details The annotation dataframe \code{ann} contains gene symbols and/or GO
#'   terms associated with the genes in the \link[DESeq2]{DESeqDataSet} object.
#'   Each row in the dataframe represents a gene. The first column contains the
#'   identifiers which correspond to the rownames of the
#'   \link[DESeq2]{DESeqDataSet} object.
#'
#'   The first column must contain the gene identifiers which correspond to the
#'   row names of the \link[DESeq2]{DESeqDataSet} object. There doesn't need to
#'   be an entry for every gene identifier in the \link[DESeq2]{DESeqDataSet}
#'   but, at least one identifier should be represented.
#'
#'   There then follows one or two further columns depending on whether the
#'   annotations consist of only gene symbols, only GO terms or both gene
#'   symbols and GO terms. Where there are multiple GO terms for a gene, the
#'   individual terms should be semi-colon (;) delimited e.g.
#'   \emph{GO:0004415;GO:0004415;GO:0004415;}.
#'
#'   \code{addAnnotationsToDataSet} returns a \link[DESeq2]{DESeqDataSet} with
#'   annotations stored in the object metadata. Gene symbols will be accessible
#'   using \code{mcols(dds)$symbol} while GO terms will be accessible using
#'   \code{metadata(dds)$go} where dds represents a \link[DESeq2]{DESeqDataSet}
#'   object.
#'
#' @family annotation functions
#'
#' @param dds A \link[DESeq2]{DESeqDataSet} object.
#' @param ann A \link{data.frame} containing gene symbols and/or delimited GO
#'   terms.
#'
#' @return A \link[DESeq2]{DESeqDataSet} object with annotations.
#'
#' @importFrom methods is
#' @importFrom GenomicRanges mcols 'mcols<-'
#' @importFrom DESeq2 DESeqDataSet
#'
#' @export
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(n=3)
#' good_ann <- data.frame( genes=rownames(dds),
#'                         gene_names=letters[1:3],
#'                         go=c("GO:0000001", "GO:0000001;GO:0000002", ""))
#' good_ann[] <- lapply(good_ann, as.character)
#' addAnnotationsToDataSet(dds, good_ann)

addAnnotationsToDataSet <- function(dds, ann)
{
  stopifnot(is(dds, 'DESeqDataSet'))
  stopifnot(is.data.frame(ann))

  if (ncol(ann) < 2) stop("Could not import annotation: too few columns in annotation.")
  if (ncol(ann) > 3) stop("Could not import annotation: too many columns in annotation.")

  col2 <- sum(grepl("GO:", ann[,2]))
  if (col2 > 0) {
    dds@metadata$go <- getGOlist(ann, 2)
  } else {
    mcols(dds)$symbol <- getGeneSymbols(dds, ann, 2)
  }

  if (ncol(ann) ==  3) {
    col3 <-  sum(grepl("GO:", ann[,3]))
    if (col2 > 0 && col3 > 0) stop("Could not import annotation: GO terms were found in multiple columns.")
    if (col2 == 0 && col3 == 0) stop("Could not import annotation: no GO terms identified.")

    annotation <- list()
    if (col2 > 0)
    {
      dds@metadata$go <- getGOlist(ann, 2)
      mcols(dds)$symbol <- getGeneSymbols(dds, ann, 3)
    }
    else if (col3 > 0)
    {
      dds@metadata$go <- getGOlist(ann,3)
      mcols(dds)$symbol <- getGeneSymbols(dds, ann, 2)
    }
  }

  return(dds)
}


#' @title Transform Delimited GO Terms To List
#'
#' @description \code{getGOlist} transforms a semi-colon delimited \link{vector}
#'   of GO terms and their corresponding gene identifiers from a
#'   \link{data.frame} into a \link{list}.
#'
#' @details The annotation dataframe \code{ann} contains gene symbols and/or GO
#'   terms associated with the genes in the \link[DESeq2]{DESeqDataSet} object.
#'   Each row in the dataframe represents a gene. The first column contains the
#'   identifiers which correspond to the rownames of the
#'   \link[DESeq2]{DESeqDataSet} object.
#'
#'   There then follows one or two further columns depending on whether the
#'   annotations consist of only gene symbols, only GO terms or both gene
#'   symbols and GO terms.  For more information see
#'   \link[deago]{annotateDataset}.
#'
#'   \code{getGOlist} requires an annotation dataframe \code{ann} and an integer
#'   value \code{col} which represents the index of the column containing the GO
#'   terms.  For example, if the GO terms were in column 2:
#'   \code{getGOlist(ann,2)}.
#'
#'   Where there are multiple GO terms for a gene, the individual terms should
#'   be semi-colon (;) delimited e.g. \emph{GO:0004415;GO:0004415;GO:0004415;}.
#'
#'   \code{getGOlist} returns a list where the names are the gene idenifiers and
#'   the values are the associated GO terms.  This format is required for
#'   \link{topGO} GO term enrichment analysis.
#'
#' @family annotation functions
#'
#' @param ann A \link{data.frame} containing gene symbols and/or delimited GO
#'   terms.
#' @param col An \link{integer} representing the index of the column containing
#'   GO terms.
#'
#' @return A \link{list} of GO terms where names are the associated gene
#'   identifiers.
#'
#' @export
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(n=3)
#' good_ann <- data.frame( genes=rownames(dds),
#'                         gene_names=letters[1:3],
#'                         go=c("GO:0000001", "GO:0000001;GO:0000002", ""))
#' good_ann[] <- lapply(good_ann, as.character)
#' getGOlist(good_ann, 3)


getGOlist <- function( ann, col )
{
  stopifnot(is.data.frame(ann))
  stopifnot(is.numeric(col) && col <= ncol(ann))

  if (sum(grepl("GO:", ann[, col])) == 0) stop("Could not import annotation: no GO terms identified.")
  go.tmp <- ann[, col]

  names(go.tmp) <- gsub(" ", "", ann[,1])
  go <- lapply(go.tmp, function(x) gsub(" ", "", strsplit(x, split=";")[[1]]))

  if (length(go) < 2) stop("Could not import annotation: GO list not created.")

  return(go)
}


#' @title Extract Gene Symbols
#'
#' @description \code{getGeneSymbols} extracts a vector of gene symbols from a
#'   dataframe.
#'
#' @details \code{getGeneSymbols} requires an annotation dataframe \code{ann}
#'   and an integer value \code{col} which represents the index of the column
#'   containing the gene symbols.  For example, if the gene symbols were in
#'   column 2: \code{getGeneSymbols(ann,2)}.
#'
#'   At least one of the gene identifiers in the row names of the
#'   \link[DESeq2]{DESeqDataSet} object must be present in the first column of
#'   the annotation dataframe \code{ann}.
#'
#'   Gene symbols are extracted from the dataframe and any missing identifiers
#'   are given a gene symbol of "unknown". There is a final check to ensure that
#'   the gene symbol vector is the same length as the row names vector from the
#'   \link[DESeq2]{DESeqDataSet} object.
#'
#'   \code{getGeneSymbols} returns a vector of gene symbols associated with the
#'   gene identifiers of the \link[DESeq2]{DESeqDataSet} object.
#'
#' @family annotation functions
#'
#' @param dds A \link[DESeq2]{DESeqDataSet} object.
#' @param ann A \link{data.frame} containing gene symbols and/or delimited GO
#'   terms.
#' @param col An \link{integer} representing the index of the column containing
#'   gene symbols.
#'
#' @return A vector of gene symbols corresponding to the rownames of the \link[DESeq2]{DESeqDataSet} object.
#'
#' @importFrom methods is
#' @importFrom DESeq2 DESeqDataSet
#'
#' @export
#'
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(n=3)
#' good_ann <- data.frame( genes=rownames(dds),
#'                         gene_names=letters[1:3],
#'                         go=c("GO:0000001", "GO:0000001;GO:0000002", ""))
#' good_ann[] <- lapply(good_ann, as.character)
#' getGeneSymbols(dds, good_ann, 2)

getGeneSymbols <- function(dds, ann, col)
{
  stopifnot(is(dds, 'DESeqDataSet'))
  stopifnot(is.data.frame(ann))
  stopifnot(is.numeric(col) && col <= ncol(ann))

  symbols <- ann[match(rownames(dds), ann[, 1]), col ]

  if (sum(!is.na(symbols)) < 1) stop("Could not import annotation: no dataset identifiers found.")

  symbols[is.na(symbols) ] <- "unknown"
  symbols[which(symbols == "") ] <- "unknown"

  if (length(symbols) != length(rownames(dds))) stop("Could not import annotation: gene symbol list length different from feature list.")

  return(symbols)
}
