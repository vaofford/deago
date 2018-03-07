#' @title Format count matrix from DESeqDataSet
#'
#' @description Format counts from DESeqDataSet for read and null counts plots
#'
#' @param dds DESeqDataSet
#'
#' @import DESeq2
#'
#' @export
#'
formatReadCounts <- function (dds)
{
  raw_counts <- counts(dds)
  normalized_counts <- counts(dds, normalized=TRUE)

  total_raw_counts <- transform(apply(raw_counts, 2, sum))
  total_normalized_counts <- transform(apply(normalized_counts, 2, sum))

  null_raw_counts <- transform(apply(raw_counts, 2, function(x) (sum(x==0)/sum(x))*100))
  null_normalized_counts <- transform(apply(normalized_counts, 2, function(x) (sum(x==0)/sum(x))*100))

  replicate_levels <- colnames(raw_counts)[order(dds$condition, colnames(raw_counts))]
  counts_df <- data.frame( "condition" = dds$condition,
                           "replicates" = factor(colnames(raw_counts), levels=replicate_levels),
                           "total_raw_counts" = total_raw_counts[,1],
                           "total_normalized_counts" = total_normalized_counts[,1],
                           "null_raw_counts" = null_raw_counts[,1],
                           "null_normalized_counts" = null_normalized_counts[,1]
                         )
  return(counts_df)
}

#' @title Plot Read Counts
#'
#' @description Generates a barplot showing the total number of reads per
#'   sample. If images subdirectory exits, will also save the plot to
#'   images/reads_per_sample.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import DESeq2
#' @import ggplot2
#' @import ggpubr
#' @importFrom scales comma
#'
#' @export

plotReadCounts <- function (dds, resultsDir)
{
  counts_df <- formatReadCounts(dds)

  ymax_raw <- roundToFactorOfTen(max(counts_df$total_raw_counts))
  rc_raw_plot <-  ggplot(counts_df, aes_string(x='replicates', y='total_raw_counts', fill='condition')) +
                  geom_bar(stat="identity") +
                  theme_deago_counts(ymax_raw) +
                  ylab("Total raw read count per sample")

  ymax_normalized <- roundToFactorOfTen(max(counts_df$total_normalized_counts))
  rc_noralized_plot <-  ggplot(counts_df, aes_string(x='replicates', y='total_normalized_counts', fill='condition')) +
                        geom_bar(stat="identity") +
                        theme_deago_counts(ymax_normalized) +
                        ylab("Total normalized read count per sample")

  rc_plot <- ggarrange( rc_raw_plot, rc_noralized_plot, ncol = 1, nrow = 2)

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(rc_plot, resultsDir, "reads_per_sample.png")
  }

  return(rc_plot)
}

#' @title Plot Null Counts
#'
#' @description Generates a barplot showing the number of genes
#'   which have no reads mapping to them (null counts)
#'   images/null_counts_per_sample.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import DESeq2
#' @import ggplot2
#' @importFrom scales comma
#'
#' @export

plotNullCounts <- function (dds, resultsDir)
{
  counts_df <- formatReadCounts(dds)

  ymax_raw <- roundToFactorOfTen(max(counts_df$null_raw_counts))
  nc_raw_plot <-  ggplot(counts_df, aes_string(x='replicates', y='null_raw_counts', fill='condition')) +
    geom_bar(stat="identity") +
    theme_deago_counts(ymax_raw) +
    ylab("Null raw counts per sample (%)")

  ymax_normalized <- roundToFactorOfTen(max(counts_df$null_normalized_counts))
  nc_noralized_plot <-  ggplot(counts_df, aes_string(x='replicates', y='null_normalized_counts', fill='condition')) +
    geom_bar(stat="identity") +
    theme_deago_counts(ymax_normalized) +
    ylab("Null normalized counts per sample (%)")

  nc_plot <- ggarrange( nc_raw_plot, nc_noralized_plot, ncol = 1, nrow = 2)

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(nc_plot, resultsDir, "null_counts_per_sample.png")
  }

  return(nc_plot)
}

#' @title Plot Sample Distances
#'
#' @description Generates a heatmap of sample-to-sample distances with cluster
#'   dendrogram. If images subdirectory exits, will also save the plot to
#'   images/sample_distances.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import DESeq2
#' @import SummarizedExperiment
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom stats dist hclust
#'
#' @export

plotSampleDistances <- function(dds, resultsDir)
{
  rld <- rlog(dds, blind=FALSE)

  sample_dists <- dist(t(assay(rld)))
  sample_order <- hclust(sample_dists)$order
  sample_levels <- colnames(rld)[sample_order]
  sample_dist_matrix <- melt(as.matrix(sample_dists))
  sample_dist_matrix$X1 <- factor(sample_dist_matrix$X1, levels=sample_levels)
  sample_dist_matrix$X2 <- factor(sample_dist_matrix$X2, levels=sample_levels)

  sd_plot <-  ggplot(data=sample_dist_matrix, aes(x=sample_dist_matrix$X1,y=sample_dist_matrix$X2)) +
              geom_tile(aes(fill=sample_dist_matrix$value)) +
              theme_deago_sample_distances()

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(sd_plot, resultsDir, "sample_distances.png", 600, 600)
  }

  return(sd_plot)
}

#' @title Calculate principal components
#'
#' @description Calcultate principal component values
#'
#' @param dds DESeqDataSet
#'
#' @import DESeq2
#' @importFrom stats prcomp
#' @importFrom genefilter rowVars
#'
#' @export

getPrincipalComponents <- function(dds) {
  rld <- rlog(dds, blind=FALSE)
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(assay(rld)[select,]))

  pc_coords <- as.data.frame(pca$x)
  pc_coords$condition <- dds$condition

  percent_var <- data.frame( pc=colnames(pca$x), proportion= pca$sdev^2 / sum( pca$sdev^2 ) )
  percent_var$cumulative <- cumsum(percent_var$proportion)

  pc_list <- list(pc_coords=pc_coords, percent_var=percent_var)
  return(pc_list)
}


#' @title PCA Plot
#'
#' @description Generates a PCA plot. If images subdirectory exits, will also
#'   save the plot to images/pca.png.
#'
#' @param pc_list List containing PC coordinates and PC summary from
#'   getPrincipalComponents(dds)
#' @param resultsDir character: path to timestamped results directory
#'
#' @import DESeq2
#' @import ggplot2
#' @import ggrepel
#'
#' @export

pcaPlot <- function (pc_list, resultsDir)
{
  data <- pc_list$pc_coords
  percent_var <- round(100*pc_list$percent_var$proportion)

  xmax = 10 * ceiling(max(data$PC1) / 10) + 10
  xmin = 10 * floor(min(data$PC1) / 10) - 10
  ymax = 10 * ceiling(max(data$PC2) / 10) + 10
  ymin = 10 * floor(min(data$PC2) / 10) - 10

  repel_label <- data$condition

  pca_plot <- ggplot(data=data, aes(x=data$PC1, y=data$PC2, color=data$condition)) +
              geom_point(size=6) +
              theme_deago_pca(percent_var[1], percent_var[2], c(xmin, xmax), c(ymin, ymax), repel_label)

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(pca_plot, resultsDir, "pca.png", 1000, 1000)
  }

  return(pca_plot)
}

#' @title PCA Scree Plot
#'
#' @description Generates a PCA scree plot. If images subdirectory exits, will also
#'   save the plot to images/pca_scree.png.
#'
#' @param pc_list List containing PC coordinates and PC summary from
#'   getPrincipalComponents(dds)
#' @param resultsDir character: path to timestamped results directory
#'
#' @import DESeq2
#' @import ggplot2
#' @importFrom stats reorder
#'
#' @export

pcaScreePlot <- function (pc_list, resultsDir)
{
  percent_var <- pc_list$percent_var

  pca_scree_plot <- ggplot( percent_var, aes(x=reorder(percent_var$pc, -percent_var$proportion), y=percent_var$proportion) ) +
    geom_bar(stat="identity") +
    geom_point(aes(x = reorder(percent_var$pc, -percent_var$proportion), y = percent_var$cumulative), shape=21, colour = "black", fill = "gray") +
    geom_hline(yintercept =0.9, linetype="dashed") +
    geom_hline(yintercept =0.7, linetype="dashed", color="gray50") +
    theme_deago_pca_scree()

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(pca_scree_plot, resultsDir, "pca_scree.png")
  }

  return(pca_scree_plot)
}

#' @title PCA Summary
#'
#' @description Generates a summary table of PCA proportions and cumulative sum
#'   values
#'
#' @param pc_list List containing PC coordinates and PC summary from
#'   getPrincipalComponents(dds)
#'
#' @import DESeq2
#' @importFrom stats reorder
#'
#' @export

pcaSummary <- function (pc_list)
{
  percent_var <- pc_list$percent_var
  percent_var[,2:3] <- round(100*percent_var[,2:3], digits=2)
  colnames(percent_var) <- c('Principal Component', 'Variation Explained (%)', 'Cumulative Total (%)')

  pca_table <- datatable( percent_var,
                          rownames = FALSE,
                          options = list(
                            pageLength = 5,
                            autoWidth = TRUE,
                            scrollCollapse = TRUE,
                            digits=3,
                            columnDefs = list(list(className = 'dt-center', targets = 0:ncol(percent_var)-1))
                         )
                        )
  return(pca_table)
}


#' @title Plot Cook's Distances
#'
#' @description Generates a boxplot showing Cook's distances per sample. If
#'   images subdirectory exits, will also save the plot to
#'   images/cooks_distances.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import ggplot2
#' @importFrom SummarizedExperiment assays
#' @importFrom reshape melt
#'
#' @export

plotCooks <- function (dds, resultsDir)
{
  replicates <- dds$condition
  names(replicates) <- colnames(dds)
  cooks_distances <- melt(as.data.frame(log10(assays(dds)[["cooks"]])))
  cooks_distances <- cooks_distances[which(!is.na(cooks_distances$value)),]

  cd_plot <-  ggplot(data=cooks_distances, aes(cooks_distances$variable, cooks_distances$value, fill=factor(replicates[cooks_distances$variable]))) +
              geom_boxplot() +
              theme_deago_cooks()

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(cd_plot, resultsDir, "cooks_distances.png", 600, 600)
  }

  return(cd_plot)
}


#' @title Plot density
#'
#' @description Generates a desnity plot. If images subdirectory exits, will
#'   also save the plot to images/normalised_density.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import ggplot2
#' @importFrom DESeq2 counts
#' @importFrom utils stack
#'
#' @export

plotDensity <- function (dds,resultsDir)
{
  toplot = data.frame(counts(dds, normalized=T))
  toplot = stack(toplot, select=colnames(toplot))
  toplot <- toplot[which(is.finite(log10(toplot$values))),]

  density_ncol <- ifelse( length(colnames(dds)) > 40, 2, 1)

  density_plot <- ggplot(toplot, aes(toplot$values, colour=toplot$ind, alpha=0.5)) +
                  geom_line(aes(color=toplot$ind), stat="density", alpha=0.5) +
                  theme_deago_density() + guides(colour=guide_legend(ncol=density_ncol))

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(density_plot, resultsDir, "normalised_density.png")
  }

  return(density_plot)
}


#' @title Plot Dispersion Estimates
#'
#' @description Generates a scatter plot showing the dispersion estimates. If
#'   images subdirectory exits, will also save the plot to
#'   images/dispersion_estimates.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import ggplot2
#' @import DESeq2
#'
#' @export

plotDispersionEstimates <- function (dds, resultsDir)
{
  px = mcols(dds)$baseMean
  sel = (px>0)
  px = px[sel]
  py = mcols(dds)$dispGeneEst[sel]
  ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  df <- data.frame(x=px, y1=pmax(py, ymin), y2=dispersions(dds)[sel], y3=mcols(dds)$dispFit[sel])

  de_plot <-  ggplot(df, aes(df$x, y = df$value, color = df$variable)) +
              geom_point(aes(y=df$y1, col="gene-est" ), cex=0.5) +
              theme_deago_dispEst()

  if (!is.null(dispersions(dds))) {
    de_plot <- de_plot + geom_point(aes(y=df$y2, col="final" ), cex=0.5)
  }

  if (!is.null(mcols(dds)$dispFit)) {
    de_plot <- de_plot + geom_point(aes(y=df$y3, col="fitted" ), cex=0.5)
  }

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(de_plot, resultsDir, "dispersion_estimates.png")
  }

  return(de_plot)
}


