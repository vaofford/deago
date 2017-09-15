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
#' @importFrom scales comma
#'
#' @export

plotReadCounts <- function (dds, resultsDir)
{
  replicates <- dds$condition
  names(replicates) <- colnames(dds)

  total_reads <- transform(apply(counts(dds), 2, sum))
  total_reads$replicates <- rownames(total_reads)
  total_reads$condition <- dds$condition
  colnames(total_reads)[1]<- c("total")

  ymax <- roundToFactorOfTen(max(total_reads$total))

  rc_plot <-  ggplot(total_reads, aes_string(x='replicates', y='total', fill='condition')) +
              geom_bar(stat="identity") +
              theme_deago_read_counts(ymax)

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(rc_plot, resultsDir, "reads_per_sample.png")
  }

  return(rc_plot)
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


#' @title PCA Plot
#'
#' @description Generates a PCA plot. If images subdirectory exits, will also
#'   save the plot to images/pca.png.
#'
#' @param dds DESeqDataSet
#' @param resultsDir character: path to timestamped results directory
#'
#' @import DESeq2
#' @import ggplot2
#' @import ggrepel
#'
#' @export

pcaPlot <- function (dds, resultsDir)
{
  rld <- rlog(dds, blind=FALSE)
  data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
  percent_var <- round(100*attr(data, "percentVar"))

  xmax = 10 * ceiling(max(data$PC1) / 10) + 10
  xmin = 10 * floor(min(data$PC1) / 10) - 10
  ymax = 10 * ceiling(max(data$PC2) / 10) + 10
  ymin = 10 * floor(min(data$PC2) / 10) - 10

  repel_label <- colnames(rld)

  pca_plot <- ggplot(data=data, aes(x=data$PC1, y=data$PC2, color=data$condition)) +
              geom_point(size=6) +
              theme_deago_pca(percent_var[1], percent_var[2], c(xmin, xmax), c(ymin, ymax), repel_label)

  image_dir <- file.path(resultsDir, "images")
  if (dir.exists(image_dir)){
    plotToFile(pca_plot, resultsDir, "pca.png", 1000, 1000)
  }

  return(pca_plot)
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


