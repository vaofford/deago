#' @title Contrast MA plot
#'
#' @description \code{\link{plotContrastMA}} generates a contrast MA plot.
#'
#' @details Add some details here.
#'
#' @param contrast A DESeqResults object
#' @param resultsDir Name of timestamped results directory
#' @param alpha FDR
#' @param lfc log2FoldChange
#' @param geneSymbols gene symbols
#' @param geneLabels description
#'
#' @import DESeq2
#' @import ggplot2
#' @import ggrepel
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export

plotContrastMA <- function(contrast, resultsDir, alpha=0.05, lfc=2, geneSymbols=NULL, geneLabels=FALSE)
{
    if (!inherits(contrast, c('matrix', 'data.frame', 'DESeqResults'))) stop ("Could not plotContrastMA: data is not a matrix, dataframe or DESeqResults object.")
    if (missing(resultsDir)) stop ("Could not plotContrastMA: no results directory given.")

    if (!is.null(metadata(contrast)$alpha)) alpha <- as.numeric(metadata(contrast)$alpha)
    if (!is.numeric(alpha) | (is.numeric(alpha) & (alpha < 0 | alpha > 1))) stop ("Could not plotContrastMA: alpha must be numeric, >0 and <1.")
    if (!is.numeric(lfc)) stop ("Could not plotContrastMA: lfc is not numeric.")

    if (is.null(geneSymbols)) geneSymbols <- `if`('symbol' %in% colnames(contrast), contrast$symbol, rownames(contrast))
    if (length(geneSymbols) != nrow(contrast)) stop("Could not plotContrastMA: length of geneSymbols is not the same as num rows in contrast.")

    baseMean <- NULL
    df <- labelDEgenes(contrast, geneSymbols, lfc, alpha)
    df <- subset(df, baseMean != 0)
    df.sig.subset <- getTopGenes(df, 5)

    contrast_name <- unlist(strsplit(deparse(substitute(contrast)), "$", fixed=TRUE))[2]

    ymax = roundToBase(df$log2FoldChange,2,'up')
    ymin = roundToBase(df$log2FoldChange,2,'down')

    ma_plot <-  ggplot(df, aes_string(x = 'baseMean', y = 'log2FoldChange')) +
                  geom_point(aes_string(color = 'de'), size = 0.8) +
                  theme_deago_ma(ymin, ymax, contrast_name, df$de)

    if (isTRUE(geneLabels)) ma_plot <- ma_plot + deago_de_labels(df.sig.subset)

    maPlotFile <- paste0(contrast_name, "_MA.png")
    image_dir <- file.path(resultsDir, "images")
    if (dir.exists(image_dir)){
      plotToFile(ma_plot, resultsDir, maPlotFile, 800, 800)
    }

    return(ma_plot)
}

#' @title Volcano plot
#' @description Generates volcano plot for each contrast
#'
#' @param contrast DESeq2 contrast
#' @param resultsDir Name of timestamped results directory
#' @param alpha FDR
#' @param lfc log2FoldChange
#' @param geneSymbols gene symbols
#' @param geneLabels description
#'
#' @import DESeq2
#' @import ggplot2
#' @import ggrepel
#' @export

plotVolcano <- function(contrast, resultsDir, alpha=0.05, lfc=2, geneSymbols=NULL, geneLabels=FALSE)
{
  if (!inherits(contrast, c('matrix', 'data.frame', 'DESeqResults'))) stop ("Could not plotVolcano: data is not a matrix, dataframe or DESeqResults object.")
  if (missing(resultsDir)) stop ("Could not plotVolcano: no results directory given.")

  if (!is.null(metadata(contrast)$alpha)) alpha <- as.numeric(metadata(contrast)$alpha)
  if (!is.numeric(alpha) | (is.numeric(alpha) & (alpha < 0 | alpha > 1))) stop ("Could not plotVolcano: alpha must be numeric, >0 and <1.")
  if (!is.numeric(lfc)) stop ("Could not plotVolcano: lfc is not numeric.")

  if (is.null(geneSymbols)) geneSymbols <- `if`('symbol' %in% colnames(contrast), contrast$symbol, rownames(contrast))
  if (length(geneSymbols) != nrow(contrast)) stop("Could not plotVolcano: length of geneSymbols is not the same as num rows in contrast.")

  df <- labelDEgenes(contrast, geneSymbols, lfc, alpha)
  df$negLogPval <- -log10(df$padj)
  df <- na.omit(df)
  df.sig.subset <- getTopGenes(df, 5)

  contrast_name <- unlist(strsplit(deparse(substitute(contrast)), "$", fixed=TRUE))[2]

  xmax = roundToBase(df$log2FoldChange,5,'up')
  xmin = roundToBase(df$log2FoldChange,5,'down')

  #log2FoldChange <- negLogPval <- NULL # To satisfy checks - should probably consider alternative solutions
  vol_plot <- ggplot(data=df, aes_string(x='log2FoldChange', y='negLogPval')) +
                geom_point(aes_string(color = 'de'), size = 0.8) +
                theme_deago_volcano(xmin, xmax, contrast_name, lfc, df$de)

  if (isTRUE(geneLabels)) vol_plot <- vol_plot + deago_de_labels(df.sig.subset)

  volPlotFile <- paste0(contrast_name, "_volcano.png")
  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(vol_plot, resultsDir, volPlotFile, 800, 800)
  }

  return(vol_plot)
}

#' @title Venn diagram
#' @description Plot venn diagram from venn count matrix
#'
#' @param countMatrix Count matrix preppared by vennCounts
#' @param resultsDir Name of timestamped results directory
#'
#' @importFrom limma vennDiagram
#' @importFrom grDevices png dev.off
#' @export

plotVenn <- function(countMatrix, resultsDir)
{
  circle_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")

  vennDiagram(countMatrix,
              show.include = FALSE,
              cex=0.9,
              circle.col=circle_colors[c(1:ncol(countMatrix-1))])

  image_dir <- file.path(resultsDir, "images")
  vennFile = file.path(image_dir, "DE_venn.png")
  if (dir.exists(image_dir)){
    png(vennFile,600,600)
    vennDiagram(countMatrix,
                show.include = FALSE,
                cex=0.9,
                circle.col=circle_colors[c(1:ncol(countMatrix-1))])
    dev.off()
  }
}


#' @title Generate count matrix for venn diagram
#' @description Generate count matrix for venn diagram for DE genes
#'
#' @param contrastList List of contrasts (pairwise comparisons)
#'
#' @importFrom limma vennCounts
#' @importFrom stats na.omit
#' @export

getVennCounts <- function ( contrastList ) {
  if ( length( names( contrastList ) ) > 5 ) return( print( "A venn diagram could not be produced as there are more than 5 contrasts." ) )
  if ( length( names( contrastList ) ) < 2 ) return( print( "A venn diagram could not be produced as there are not enough contrasts." ) )

  alpha <- ifelse( is.null( metadata(contrastList[[1]])$alpha ), 0.05, metadata(contrastList[[1]])$alpha)
  alpha <- as.numeric(alpha)

  sig_gene_df <- list()
  for ( contrast in names(contrastList) )
  {
    cont <- as.data.frame(contrastList[[contrast]])
    cont_filt <- na.omit( cont )
    cont_filt <- cont_filt[(cont_filt$log2FoldChange > 2 | cont_filt$log2FoldChange < -2) & cont_filt$padj < alpha,]
    sig_genes <- rownames( cont_filt )
    cont_df <- data.frame( row.names=sig_genes, present=rep( 1, length( sig_genes ) ) )
    sig_gene_df[[contrast]] <- cont_df
  }

  all_sig_genes <- merge(sig_gene_df[1],sig_gene_df[2],by="row.names",all=T, row.names='Row.names')

  if (length( sig_genes ) > 2 ) {
    for ( i in 3:length( sig_genes ) ) {
      tmp.merge <- merge(all_sig_genes,sig_gene_df[i],by.x="Row.names", by.y="row.names",all=T, row.names="Row.names")
      all_sig_genes <- tmp.merge
    }
  }

  all_sig_genes[is.na(all_sig_genes)] <- 0
  colnames(all_sig_genes) <- c("GeneID", names(contrastList))
  all_counts <- vennCounts(all_sig_genes[, -which( colnames(all_sig_genes) == "GeneID" )])

  return(all_counts)
}

#' @title Subset top DE genes
#' @description Helper function to extract n DE genes with the greatest
#'   Log2FoldChange and n genes with the lowest Log2FoldChange from a
#'   DESeqResults dataframe.
#'
#' @param df Dataframe from a \link{DESeqResults} object
#' @param n  Number of genes to label in each direction
#'
#' @importFrom utils head tail
#'
#' @export
#'
getTopGenes <- function(df,n)
{
  df_sig <- df[which(df$de != 0), ]
  df_sig <- df_sig[order(df_sig$log2FoldChange, decreasing = TRUE), ]

  if (nrow(df_sig) > (2*n) )
  {
    df_sig_subset <- rbind(head(df_sig,n),tail(df_sig,n))
  } else {
    df_sig_subset <- df_sig
  }

  return(df_sig_subset)
}

#' @title Label DE genes
#' @description Helper function to label DE genes from a DESeqResults dataframe
#'   with a given Log2FoldChange and padj cutoff. Adds a label column (de) where
#'   values are non-significant (0), up-regulated (1) or down-regulated (2).
#'
#'
#' @param contrast A \link{DESeqResults} object
#' @param geneSymbols details
#' @param lfc Log2FoldChange threshold value
#' @param alpha FDR threshold
#'
#' @export
labelDEgenes <- function(contrast, geneSymbols=NULL, lfc=2, alpha=0.05)
{
  if (!inherits(contrast, c('matrix', 'data.frame', 'DESeqResults'))) stop ("Could not plotContrastMA: data is not a matrix, dataframe or DESeqResults object.")
  if (is.null(geneSymbols)) geneSymbols <- `if`('symbol' %in% colnames(contrast), contrast$symbol, rownames(contrast))
  if (length(geneSymbols) != nrow(contrast)) stop("Could not plotContrastMA: length of geneSymbols is not the same as num rows in contrast.")

  df <- data.frame(genes=rownames(contrast), baseMean=contrast$baseMean, log2FoldChange=contrast$log2FoldChange, padj=contrast$padj, de=rep(0, nrow(contrast)), geneSymbols=geneSymbols)
  df$de[which(df$log2FoldChange >= lfc & df$padj< alpha)] <- 1
  df$de[which(df$log2FoldChange <= -lfc & df$padj < alpha)] <- 2
  df$de <- as.factor(df$de)

  return(df)
}
