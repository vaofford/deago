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

  gene_levels <- c("all", "up", "down")

  go_table_list <- list()
  for (contrast in names(contrasts))
  {
    for (go_level in go_levels)
    {
      for (gene_level in gene_levels) {
        
        if ( gene_level != "all" ) {
          go_label <- paste(contrast, go_level, gene_level, sep="_" )
        } else {
          go_label <- paste0(contrast, "_", go_level)
        }
        
        go_data <- prepareGOdata(dds, contrasts[[contrast]], go_level, gene_level)
        
        if ( class(go_data) == 'topGOdata' ) {
          go_table <- topGOanalysis(go_data)
    
          go_table$identifiers <- getGOidentifiers(contrasts[[contrast]], go_data, go_table, gene_level)
    
          if ("symbol" %in% names(contrasts[[contrast]])) {
            go_table$symbol <- getGOsymbols(contrasts[[contrast]], go_data, go_table, gene_level)
          }
    
          go_table_list[[go_label]] <- go_table
        } else {
          go_table_list[[go_label]] <- go_data
        }
      }
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
#' @param gene_level all, up or down [all]
#'
#' @import topGO
#' @importFrom methods new
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export

prepareGOdata <- function(dds, contrast, go_level, gene_level='all')
{
  alpha <- ifelse(!is.null(metadata(contrast)$alpha), 0.05, metadata(contrast)$alpha)
  alpha <- as.numeric(alpha)
  if (!is.numeric(alpha) | (is.numeric(alpha) & (alpha < 0 | alpha > 1))) stop ("Could not prepareGoData: alpha must be numeric, >0 and <1.")

  all_genes <- rownames(contrast)
  contrast_filt <- as.data.frame(contrast[which(contrast$padj < alpha & !is.na(contrast$padj)),])
  
  if ( gene_level == "up") {
    contrast_filt_up <- contrast_filt[ which( contrast_filt$log2FoldChange > 0 ) , ]
    sig_genes <- rownames(contrast_filt_up)
  } else if (gene_level == "down") {
    contrast_filt_down <- contrast_filt[ which( contrast_filt$log2FoldChange < 0 ) , ]
    sig_genes <- rownames(contrast_filt_down)
  } else { # Assume all 
    sig_genes <- rownames(contrast_filt)
  }

  relevant_genes <- factor(as.integer(all_genes %in% sig_genes))
  names(relevant_genes) <- all_genes

  groupGOTerms() # Load the GOTerm environments (as of 3.5)
  
  if ( length(sig_genes) > 0 ) {
    go_data <-new("topGOdata",
                 ontology = go_level,
                 allGenes = relevant_genes,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = metadata(dds)$go
    )
  } else {
    return( "A GO enrichment analysis could not be performed as there are no significantly differentially expressed genes for this contrast." )
  }

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

#' @title Get gene identifiers for GO terms
#' @description Get gene identifiers for GO terms
#' 
#' @param contrast DESeqResults object
#' @param GOdata topGO object
#' @param GOtable topGO table
#' @param gene_level all, up or down [all]
#'
#' @import topGO
#' @export

getGOidentifiers <- function(contrast, GOdata, GOtable, gene_level='all')
{
  goGeneList <- genesInTerm(GOdata, GOtable$GO.ID)
  filteredGoGeneList <- filterGeneIdentifiers(contrast, goGeneList, gene_level)

  identifiersInTerms <- vector()
  for (i in 1:length(GOtable$GO.ID))
  {
    goTerm <- GOtable$GO.ID[i]
    genesInTerm <- filteredGoGeneList[goTerm][[1]]
    identifiersInTerms[i] <- paste(genesInTerm, collapse=', ')
  }

  return(identifiersInTerms)
}

#' @title Get gene symbols for GO terms
#' @description Get gene symbols for GO terms
#'
#' @param contrast DESeqResults object
#' @param GOdata topGO object
#' @param GOtable topGO table
#' @param gene_level all, up or down [all]
#'
#' @import topGO
#' @export

getGOsymbols <- function(contrast, GOdata, GOtable, gene_level='all')
{
  goGeneList <- genesInTerm(GOdata, GOtable$GO.ID)
  filteredGoGeneList <- filterGeneIdentifiers(contrast, goGeneList, gene_level)

  symbolsInTerms <- vector()
  for (i in 1:length(GOtable$GO.ID))
  {
    goTerm <- GOtable$GO.ID[i]
    genesInTerm <- filteredGoGeneList[goTerm][[1]]
    geneIndex <- rownames(contrast) %in% genesInTerm
    symbolsInTerm <- contrast$symbol[geneIndex]
    symbolsInTerm <- sort( unique(symbolsInTerm), method='radix' )
    symbolsInTerms[i] <- paste(symbolsInTerm, collapse=', ')
  }

  return(symbolsInTerms)
}

#' @title Filter DE gene identifiers in GO term
#' @description Filter DE gene identifiers for GO terms
#' 
#' @param contrast DESeqResults object
#' @param goGeneList vector of gene identifiers associated with GO term
#' @param gene_level all, up or down [all]
#'
#' @import topGO
#' @export

filterGeneIdentifiers <- function ( contrast, goGeneList, gene_level='all' )
{ 
  if ( gene_level == 'up' ) {
    genes <- rownames(contrast)[ which( contrast$padj < 0.05 & contrast$log2FoldChange > 0 ) ]
  } else if ( gene_level == 'down' ) {
    genes<- rownames(contrast)[ which( contrast$padj < 0.05 & contrast$log2FoldChange < 0 ) ]
  } else {
    genes <- rownames(contrast)[ which(contrast$padj < 0.05) ]
  }
  
  go.genes <- list()
  for ( go.id in names(goGeneList) ) {
    go.genes[[go.id]] <- intersect( genes, unlist(goGeneList[go.id]) )
  }
  
  return(go.genes)
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
  if ( !is.data.frame(go_table) ) {
    message("A GO enrichment analysis could not be performed as there are no significantly differentially expressed genes for this contrast.")
  } else {
    condensed_go_table <- go_table[, c('GO.ID','Term','Significant','Expected','weight01Fisher')]
    if("symbol" %in% colnames(go_table))
    {
      condensed_go_table$symbol <- go_table$symbol
    }
    
    condensed_go_table$GO.ID <- paste0("<a href='http://amigo.geneontology.org/amigo/term/",
                                       condensed_go_table$GO.ID, "'target='_blank'>",
                                       condensed_go_table$GO.ID,"</a>")
  
    go_dt <- datatable(condensed_go_table,
                       filter = 'top',
                       escape=FALSE,
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
}
