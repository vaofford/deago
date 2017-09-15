#' @title Extract contrasts from DESeq2 results
#'
#' @description Extract contrasts from DESeq2 results and write to file
#'
#' @details
#' The details
#'
#' @param dds A \link[DESeq2]{DESeqDataSet} object.
#' @param parameters A \link{list} containing key/value pairs which define the
#'   parameters for the analysis (see \link[deago]{importConfig} for more
#'   information).
#'
#' @importFrom methods is
#' @import DESeq2
#' @importFrom utils combn write.table
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export

getContrasts <- function(dds, parameters)
{
  stopifnot(is(dds, 'DESeqDataSet'))
  stopifnot(is.list(parameters))

  alpha <- ifelse(is.null(parameters$alpha), 0.05, parameters$alpha)

  # Generate contrasts
  contrast_matrix <- combn(levels(dds$condition), 2)
  if(ncol(contrast_matrix) > 1)
  {
      contrasts <- apply(contrast_matrix[c(2, 1), ] , 2 , paste , collapse = "_vs_")
  } else {
      contrasts <- paste(contrast_matrix[2,], contrast_matrix[1,], sep='_vs_')
  }

  results <- list()
  for (i in 1:ncol(contrast_matrix))
  {
    # Extract contrast from results
    results[[contrasts[i]]] <- results(dds, alpha=alpha, contrast=c("condition", contrast_matrix[2, i], contrast_matrix[1, i]))

    if ( "symbol" %in% names( mcols( dds ) ) ) {
      results[[contrasts[i]]]$symbol <- mcols(dds)$symbol
    }

    if ( "go" %in% names(metadata(dds))) {
      metadata(results[[contrasts[i]]])$go <- metadata(dds)$go
    }
  }

  return(results)
}

#' @title Extract contrasts from DESeq2 results
#'
#' @description Extract contrasts from DESeq2 results and write to file
#'
#' @details
#' The details
#'
#' @param dds A \link[DESeq2]{DESeqDataSet} object.
#' @param contrasts A \link{list}
#' @param resultsDir A \link{character} string giving the path to timestamped
#'   results directory
#'
#' @importFrom methods is
#' @import DESeq2
#' @importFrom utils combn write.table
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export

writeContrasts <- function(dds, contrasts, resultsDir)
{
  for (i in names(contrasts))
  {
    alpha <- ifelse(is.null(metadata(contrasts[[i]])$alpha), 0.05, metadata(contrasts[[i]])$alpha)
    results_file <- paste0(i,"_q",alpha,".txt")
    results_table <- prepareContrast(dds, contrasts[[i]])
    write.table(results_table, file=file.path(resultsDir, results_file), quote=FALSE, sep="\t", row.names=FALSE)
  }
}

#' @title Prepares contrasts table for writing
#' @description Takes results object and prepares dataframe for writing
#'
#' @param dds DESeqDataSet object containing DESeq2 analysis results
#' @param contrast DESeqResults object
#'
#' @import DESeq2
#' @importFrom methods as
#' @export

prepareContrast <- function(dds, contrast)
{
  normalised_counts <- counts(dds, normalized=T)
  contrast <- cbind(normalised_counts, as( contrast, "matrix" ) ) # Need to transform to matrix or throws error

  # Add gene symbol if available
  if ( "symbol" %in% names( mcols( dds ) ) )
  {
    contrast <- cbind( mcols(dds)$symbol, contrast )
    colnames(contrast)[1] <- "symbol"
  }

  # Horrible hack to stop shift of column names when keeping column headers
  geneIDs <- rownames(contrast)
  contrast <- cbind(geneIDs, contrast)
  colnames(contrast)[1] <- "geneID"

  return(contrast)
}


#' @title Differential expression summary
#' @description Summary table of differentially expressed gene number per contrast
#'
#' @param contrastList list of contrast results objects
#' @param parameters list of parameters
#'
#' @export

contrastSummary <- function(contrastList, parameters)
{
  alpha <- ifelse(is.null(parameters$alpha), 0.05, parameters$alpha)

  contrast_table <- data.frame(row.names=names(contrastList))

  for ( contrast in names(contrastList) ) {
    cont <- contrastList[[contrast]]
    up <- length( which( cont$log2FoldChange > 2 & cont$padj < abs(alpha) ) )
    down <- length( which( cont$log2FoldChange < -2 & cont$padj < abs(alpha) ) )
    total <- up + down
    contrast_table[contrast,1] <- up
    contrast_table[contrast,2] <- down
    contrast_table[contrast,3] <- total
  }

  names(contrast_table) <- c("up-regulated", "down-regulated", "total")

  return(contrast_table)
}


#' @title Prepare contrast table
#' @description Prepare a filtered table of contrast results containing DE genes only
#'
#' @param contrast DESeqResults object
#'
#' @import DESeq2
#' @import DT
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export
#'
prepareContrastTable <- function(contrast)
{
  contrast_filt <- contrast[which(contrast$padj < 0.01 & !is.na(contrast$padj) & (contrast$log2FoldChange <= -2 | contrast$log2FoldChange >= 2) ),]

  if (nrow(contrast_filt) < 1)
  {
    contrast_name <- unlist(strsplit(deparse(substitute(contrast)), "$", fixed=TRUE))[2]
    message(paste0("Could not generate contrast table for ", contrast_name, " as no differentially expressed genes were found (q<0.01)."))
  } else {
    contrast_filt <- as.data.frame(contrast_filt)[,c('baseMean','log2FoldChange','padj')] %>% format(digits=3)
    contrast_filt[,1:2] <- apply(contrast_filt[,1:2], 2, as.numeric)

    if ( "symbol" %in% names( contrast ) )
    {
      symbols <-  contrast$symbol
      symbols_filt <- symbols[which(contrast$padj < 0.01 & !is.na(contrast$padj) & (contrast$log2FoldChange <= -2 | contrast$log2FoldChange >= 2) )]
      contrast_filt <- cbind(symbols_filt, contrast_filt )
      colnames(contrast_filt)[1] <- "symbol"
    }

    contrast_filt <- contrast_filt[order(contrast_filt$padj),]

    if (substr(rownames(contrast)[1],1,3) == "ENS") {
      rownames(contrast_filt) <- paste0("<a href='http://www.ensembl.org/id/",
                                         rownames(contrast_filt), "'target='_blank'>",
                                         rownames(contrast_filt),"</a>")
    }

    contrast_table <- datatable(contrast_filt,
                        filter = 'top',
                        escape=FALSE,
                        options = list(
                          pageLength = 10,
                          autoWidth = TRUE,
                          scrollX = TRUE,
                          scrollCollapse = TRUE,
                          digits=3,
                          columnDefs = list(list(className = 'dt-center', targets = 0:ncol(contrast_filt)-1))
                      )
    )
    return(contrast_table)
  }
}
