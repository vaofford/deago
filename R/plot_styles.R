theme_deago <- function()
{
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    axis.title=element_text(size=14, face="bold", margin=margin(0, 15, 0, 0)),
    plot.title=element_text(size=16, face="bold", hjust=0.5)
  )
}

theme_deago_legend <- function()
{
  theme(
    legend.title=element_text(size=14, face="bold", vjust=2),
    legend.text=element_text(size=12)
  )
}

theme_deago_density_legend <- function()
{
  theme_deago_legend() +
  theme(
    legend.justification=c(1,0),
    legend.position=c(1,0),
    legend.background = element_rect(fill="gray90", size=.5, linetype="solid", color = "black")
  )
}

theme_deago_de_legend <- function()
{
  theme_deago() +
  theme(
    legend.position="top",
    legend.title=element_blank(),
    legend.text=element_text(size = 14, face = "bold")
  )
}

scale_color_de <- function(sig)
{
  ns.lab <- paste("ns = ", length(sig))
  up.lab <- paste("up = ", sum(sig==1))
  down.lab <- paste("down = ", sum(sig==2))

  de_palette <- c('1'='#D55E00', '2'='#009E73','0'='lightgray')
  structure(list(
    scale_color_manual(values=de_palette, labels = c(ns.lab, up.lab, down.lab)),
    guides(colour = guide_legend(override.aes = list(size = 3),
                                 keywidth = 0.5, keyheight = 0.1, default.unit = "inch"))
  ))
}

deago_de_labels <- function(topGenes)
{
  structure(list(
    geom_label_repel(data = topGenes, mapping = aes_string(label = 'geneSymbols'),
                     box.padding = unit(1, "lines"),
                     size = 5, fontface = 'italic')
  ))
}

theme_deago_ma <- function(ymin, ymax, title, sig)
{
  structure(list(
    scale_x_log10(breaks = c(0.1, 10, 1000, 100000)),
    ylim(ymin,ymax),
    labs( x="mean of normalized counts",
          y="log fold change",
          title=title),
    theme_deago_de_legend(),
    scale_color_de(sig)
  ))
}

theme_deago_volcano <- function(xmin, xmax, title, lfc, sig)
{
  structure(list(
    xlim(xmin,xmax),
    labs( x="log2FoldChange",
          y="pvalue(-log10)",
          title=title),
    theme_deago_de_legend(),
    scale_color_de(sig),
    geom_vline(xintercept=lfc),
    geom_vline(xintercept=-lfc)
  ))
}

theme_deago_read_counts <- function(ymax)
{
  structure(list(
    theme_deago(),
    theme_deago_legend(),
    xlab(""),
    ylab("Total read count per sample"),
    scale_fill_discrete(name="Experimental condition"),
    theme(axis.text.x = element_text(angle = 90, hjust = 1)),
    scale_y_continuous(labels=comma, limits=c(0, ymax))
  ))
}

theme_deago_sample_distances <- function(ymax)
{
  structure(list(
    theme_deago(),
    theme_deago_legend(),
    xlab(""),
    ylab(""),
    scale_fill_gradient(low="darkorange4", high="cornsilk", name="Distance"),
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ))
}

theme_deago_pca <- function(pc1, pc2, xlim, ylim, repel_label)
{
  structure(list(
    theme_deago(),
    theme_deago_legend(),
    xlab(paste0("PC1: ", pc1, "% variance")),
    ylab(paste0("PC2: ", pc2, "% variance")),
    xlim(xlim),
    ylim(ylim),
    scale_color_discrete(name = "Experimental condition"),
    geom_label_repel(aes(label=repel_label, size=3), show.legend = FALSE, box.padding = unit(1, "lines"))
  ))
}

theme_deago_cooks <- function()
{
  structure(list(
    theme_deago(),
    theme_deago_legend(),
    xlab(""),
    ylab("Cook's distance (log10)"),
    scale_fill_discrete(name = "Experimental condition"),
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ))
}

theme_deago_density <- function()
{
  structure(list(
    theme_deago(),
    theme_deago_legend(),
    scale_x_log10(name="\nNormalized counts",
                  breaks=c(0.1,1,10,100,1000,10000,100000),
                  limits=c(0.1,100000)),
    scale_y_continuous(name="Density\n"),
    scale_colour_discrete(name="Samples"),
    geom_vline(xintercept=10, colour="grey", linetype = "dashed")
  ))
}

theme_deago_dispEst <- function()
{
  structure(list(
    theme_deago(),
    theme_deago_density_legend(),
    scale_x_log10(),
    scale_y_log10(),
    xlab("mean of normalized counts"),
    ylab("dispersion"),
    scale_color_manual(name="Legend", values=c("gene-est"="black", "final"="dodgerblue","fitted"="red"))
  ))
}
