library(deago)
library(DESeq2)

expected_parameters_list <- list (  'counts_directory'  = "test_counts",
                                    'targets_file'      = "test_targets.txt",
                                    'results_directory' = "/path/place/holder",
                                    'gene_ids'          = "Geneid",
                                    'alpha'             = 0.05,
                                    'control'           = NULL,
                                    'columns'           = 'condition',
                                    'annotation_file'   = 'tests/testthat/deago-test-annotation.tsv',
                                    'keep_images'       = 1,
                                    'qc_only'           = 0,
                                    'go_analysis'       = 1)

parameters_file <- file.path(tempdir(), 'config')
buildConfig(parameters_file, parameters=expected_parameters_list)
expected_parameters <- importConfig(parameters_file)
unlink(parameters_file)

expected_dds <- makeExampleDESeqDataSet(m=18, betaSD=2)
expected_dds$condition <- factor(c(rep(c("AI","AII","AIII"),each=3), rep(c("BI","BII","BIII"),each=3)))

target_df <- data.frame(  "condition" = as.character(expected_dds$condition),
                          "replicate" = rep(c("1", "2", "3"), times=6),
                          "filename" = paste0(colnames(expected_dds),".tsv"),
                          stringsAsFactors=FALSE)
targets_file <- tempfile("test_targets", fileext=".txt")
write.table(target_df, file=targets_file, quote=FALSE, row.names=FALSE, sep="\t")
expected_targets <- importTargets(targets_file)
unlink(targets_file)

expected_dds_de <- suppressMessages(DESeq(expected_dds))
colData(expected_dds)$condition <- factor(tolower(colData(expected_dds)$condition))

expected_annotation <- importAnnotation(expected_parameters$annotation_file)
expected_dds_de_ann <- annotateDataset(expected_dds_de, expected_parameters)
expected_gene_list <- getGeneSymbols(expected_dds_de, expected_annotation, 2)
expected_go_list <- getGOlist(expected_annotation, 3)

expected_rc_plot <- plotReadCounts(expected_dds_de, getwd())
expected_sd_plot <- plotSampleDistances(expected_dds_de, getwd())
expected_pca_plot <- pcaPlot(expected_dds_de, getwd())
expected_cooks_plot <- plotCooks(expected_dds_de, getwd())
expected_density_plot <- plotDensity(expected_dds_de, getwd())
expected_dispersion_plot <- plotDispersionEstimates(expected_dds_de, getwd())

expected_contrasts <- getContrasts(expected_dds_de, expected_parameters)
expected_contrasts_ann <- getContrasts(expected_dds_de_ann, expected_parameters)
expected_prepared_contrast <- prepareContrast(expected_dds_de, expected_contrasts[["BI_vs_AI"]])
expected_contrast_summary <- contrastSummary(expected_contrasts, list())
expected_contrast_table <- prepareContrastTable(expected_contrasts[["BI_vs_AI"]])

expected_go_tables <- runGOanalysis(expected_dds_de_ann, list("BI_vs_AI"=expected_contrasts_ann[["BI_vs_AI"]]) , c('BP'))
expected_go_data <- prepareGOdata(expected_dds_de_ann, expected_contrasts_ann[["BI_vs_AI"]], 'BP')
expected_go_table <- topGOanalysis(expected_go_data)
expected_go_symbols <- getGOsymbols(expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table)
expected_go_dt <- prepareGOtable(expected_go_tables[["BI_vs_AI_BP"]])

expected_ma_plot <- plotContrastMA(expected_contrasts$BI_vs_AI, getwd(), geneLabels=TRUE)
expected_volcano_plot <- plotVolcano(expected_contrasts$BI_vs_AI, getwd(), geneLabels=TRUE)
expected_venn_counts <- getVennCounts(expected_contrasts[c(3,12)])

gene_symbols <- paste0("test_", rownames(expected_contrasts[[3]]))
expected_labelled_genes <- labelDEgenes(expected_contrasts[[3]], geneSymbol=gene_symbols, lfc=3, alpha=0.01)
expected_top_genes <- getTopGenes(expected_labelled_genes, 10)

sessInfo <- sessionInfo()

save( expected_parameters,
      expected_dds,
      target_df,
      expected_targets,
      expected_dds_de,
      expected_annotation,
      expected_dds_de_ann,
      expected_gene_list,
      expected_go_list,
      expected_contrasts,
      expected_contrasts_ann,
      expected_prepared_contrast,
      expected_contrast_summary,
      expected_contrast_table,
      expected_go_tables,
      expected_go_data,
      expected_go_table,
      expected_go_symbols,
      expected_go_dt,
      expected_rc_plot,
      expected_sd_plot,
      expected_pca_plot,
      expected_cooks_plot,
      expected_density_plot,
      expected_dispersion_plot,
      expected_ma_plot,
      expected_volcano_plot,
      expected_venn_counts,
      expected_labelled_genes,
      expected_top_genes,
      sessInfo,
      file="tests/testthat/deago-testdata.RData")





