#' @title GO enrichment analysis
#' @description Run topGO enrichment anlaysis for contrast
#'
#' @param dds DESeq object
#' @param contrast DESeqResults object
#' @param go_level BP or MF
#'
#' @import topGO
#' @importFrom S4Vectors metadata 'metadata<-'
#'
#' @export

prepareGOtable <- function(dds, contrast, go_level)
{
  go_data <- prepareGoData(dds, contrast, go_level)

  go_table <- runGoAnalysis(go_data)

  if ("symbol" %in% names(contrast)) {
    go_table <- getGoSymbols(contrast, go_data, go_table)
  }

  go_dt <- datatable(go_table,
                     filter = 'top',
                     options = list(
                       pageLength = 10,
                       autoWidth = TRUE,
                       scrollX = TRUE,
                       scrollCollapse = TRUE,
                       digits=3,
                       rownames= FALSE,
                       columnDefs = list(list(className = 'dt-center', targets = 1:ncol(go_table)-1))
                     )
                    )

  return(go_dt)
}

#' @title Prepare GO data
#' @description Prepare GO terms for TopGO analysis
#'
#' @param dds DESeq object
#' @param contrast DESeqResults object
#' @param go_level BP or MF
#'
#' @import topGO
#' @importFrom methods new
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export

prepareGoData <- function(dds, contrast, go_level)
{
  alpha <- ifelse(!is.null(metadata(contrast)$alpha), 0.05, metadata(contrast)$alpha)
  alpha <- as.numeric(alpha)
  if (!is.numeric(alpha) | (is.numeric(alpha) & (alpha < 0 | alpha > 1))) stop ("Could not prepareGoData: alpha must be numeric, >0 and <1.")

  all_genes <- rownames(contrast)
  contrast_filt <- as.data.frame(contrast[which(contrast$padj < alpha & !is.na(contrast$padj)),])
  sig_genes <- rownames(contrast_filt)
  relevant_genes <- factor(as.integer(all_genes %in% sig_genes))
  names(relevant_genes) <- all_genes

  groupGOTerms() # Load the GOTerm environments (as of 3.5)

  go_data <-new("topGOdata",
               ontology = go_level,
               allGenes = relevant_genes,
               nodeSize = 10,
               annot = annFUN.gene2GO,
               gene2GO = metadata(dds)$go
  )

  return(go_data)
}

#' @title Prepare GO data
#' @description Prepare GO terms for TopGO analysis
#'
#' @param go_data topGO object
#'
#' @import DESeq2
#' @import topGO
#' @importFrom graph numNodes
#' @export

runGoAnalysis <- function(go_data)
{
  fisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  weight01 <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  elim <- runTest(go_data, algorithm = "elim", statistic = "fisher")

  minNodes <- numNodes(graph(go_data))
  numTableNodes <- 30

  if (minNodes < numTableNodes) {
    numTableNodes <- minNodes
  }

  go_table <- GenTable(go_data,
                       classicFisher = fisher,
                       elimFisher = elim,
                       weight01Fisher = weight01,
                       orderBy = "weight01Fisher",
                       ranksOf = "classicFisher",
                       numChar = 100,
                       topNodes = numTableNodes
  )

  go_table <- go_table[, c('GO.ID','Term','Significant','Expected','weight01Fisher')]

  return(go_table)
}


#' @title Get gene symbols for GO terms
#' @description Get gene symbols for GO terms
#'
#' @param contrast DESeqResults object
#' @param GOdata topGO object
#' @param GOtable topGO table
#'
#' @import topGO
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export

getGoSymbols <- function(contrast, GOdata, GOtable)
{
  goGeneList <- genesInTerm(GOdata, GOtable$GO.ID)

  for (i in 1:length(GOtable$GO.ID))
  {
    goTerm <- GOtable$GO.ID[i]
    genesInTerm <- goGeneList[goTerm][[1]]
    geneIndex <- rownames(contrast) %in% genesInTerm
    symbolsInTerm <- contrast$symbol[geneIndex]
    symbolsInTerm <- sort( unique(symbolsInTerm), method='radix' )
    GOtable$symbol[i] <- paste(symbolsInTerm, collapse=', ')
  }

  return(GOtable)
}
