#' @title GO enrichment analysis
#' @description Run topGO enrichment anlaysis for contrast
#'
#' @param dds DESeq object
#' @param contrasts A \link{list}
#' @param parameters A \link{list} containing key/value pairs which define the
#'   parameters for the analysis (see \link[deago]{importConfig} for more
#'   information).
#'
#' @import topGO
#' @importFrom S4Vectors metadata 'metadata<-'
#'
#' @export

runGOanalysis <- function(dds, contrasts, parameters)
{
  allowed_go_levels <- list('BP'=c('BP'), 'MF'=c('MF'), 'all'=c('BP','MF'))
  go_levels <- ifelse(is.null(parameters$go_level), 'all', parameters$go_level)

  if ( go_levels %in% names(allowed_go_levels) )
  {
    go_levels <- allowed_go_levels[[go_levels]]
  } else {
    stop ("go_level must be BP, MF or all")
  }

  go_table_list <- list()
  for (contrast in names(contrasts))
  {
    for (go_level in go_levels)
    {
      go_label <- paste0(contrast, "_", go_level)

      go_data <- prepareGOdata(dds, contrasts[[contrast]], go_level)
      go_table <- topGOanalysis(go_data)

      go_table$identifiers <- paste(genesInTerm(go_data, go_table$GO.ID), collapse=", ")

      if ("symbol" %in% names(contrasts[[contrast]])) {
        go_table <- getGOsymbols(contrasts[[contrast]], go_data, go_table)
      }

      go_table_list[[go_label]] <- go_table
    }
  }
  return(go_table_list)
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

prepareGOdata <- function(dds, contrast, go_level)
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

#' @title Run topGO analysis
#' @description Run topGO analysis
#'
#' @param go_data topGO object
#'
#' @import DESeq2
#' @import topGO
#' @importFrom graph numNodes
#' @export

topGOanalysis <- function(go_data)
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
#' @export

getGOsymbols <- function(contrast, GOdata, GOtable)
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

#' @title Write GO results table to file
#' @description Write GO results table to file
#'
#' @details
#' The details
#'
#' @param go_tables A list containing topGO results tables
#' @param resultsDir A \link{character} string giving the path to timestamped
#'   results directory
#'
#' @importFrom methods is
#' @importFrom utils write.table
#' @export

writeGOtables <- function(go_tables, resultsDir) {
  for (go_table in names(go_tables))
  {
    go_file = paste0(go_table, ".tsv")
    write.table(go_tables[[go_table]], file=file.path(resultsDir, go_file), quote=FALSE, sep="\t", row.names=FALSE)
  }
}


#' @title Prepare GO results datatable
#' @description Prepare GO results datatable
#'
#' @param go_table A data frame containing topGO results
#'
#' @import DT
#' @export

prepareGOtable <- function(go_table)
{
  condensed_go_table <- go_table[, c('GO.ID','Term','Significant','Expected','weight01Fisher')]
  if("symbol" %in% colnames(go_table))
  {
    condensed_go_table$symbol <- go_table$symbol
  }

  go_dt <- datatable(condensed_go_table,
                     filter = 'top',
                     options = list(
                       pageLength = 10,
                       autoWidth = TRUE,
                       scrollX = TRUE,
                       scrollCollapse = TRUE,
                       digits=3,
                       rownames= FALSE,
                       columnDefs = list(list(className = 'dt-center', targets = 1:ncol(condensed_go_table)-1))
                     )
  )

  return(go_dt)
}
