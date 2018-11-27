context("contrast plots")

test_that("plotting MA works", {
  expect_silent(observed_ma_plot <- plotContrastMA(expected_contrasts$BI_vs_AI, expected_parameters$results_directory, geneLabels=TRUE))
  expect_is(observed_ma_plot, 'ggplot')

  observed_ma_file <- file.path(expected_parameters$results_directory, 'images', 'BI_vs_AI_MA.png')
  expect_true(file.exists(observed_ma_file))

  expect_identical(layer_data(observed_ma_plot), layer_data(expected_ma_plot))
  expect_equal(layer_scales(observed_ma_plot), layer_scales(expected_ma_plot))
  expect_identical(observed_ma_plot$theme, expected_ma_plot$theme)

  expect_silent(observed_ma_plot <- plotContrastMA(expected_contrasts_ann$BI_vs_AI, expected_parameters$results_directory, geneLabels=TRUE))
})

test_that("plotting MA with bad parameters fails", {
  expect_error(plotContrastMA(list()), "Could not plotContrastMA: data is not a matrix, dataframe or DESeqResults object.")
  expect_error(plotContrastMA(matrix()), "Could not plotContrastMA: no results directory given.")
  
  bad_contrast <- expected_contrasts_ann$BI_vs_AI
  metadata(bad_contrast)$alpha <- 2
  expect_error(plotContrastMA(bad_contrast, expected_parameters$results_directory), "Could not plotContrastMA: alpha must be numeric, >0 and <1.")
  
  expect_error(plotContrastMA(expected_contrasts_ann$BI_vs_AI, expected_parameters$results_directory, lfc='x'), "Could not plotContrastMA: lfc is not numeric.")
  expect_error(plotContrastMA(expected_contrasts_ann$BI_vs_AI, expected_parameters$results_directory, geneSymbols=c('x')), "Could not plotContrastMA: length of geneSymbols is not the same as num rows in contrast.")
})

test_that("plotting volcano works", {
  expect_silent(observed_volcano_plot <- plotVolcano(expected_contrasts$BI_vs_AI, expected_parameters$results_directory, geneLabels=TRUE))
  expect_is(observed_volcano_plot, 'ggplot')

  observed_volcano_file <- file.path(expected_parameters$results_directory, 'images', 'BI_vs_AI_volcano.png')
  expect_true(file.exists(observed_volcano_file))

  expect_equal(layer_data(observed_volcano_plot), layer_data(expected_volcano_plot))
  expect_equal(layer_scales(observed_volcano_plot), layer_scales(expected_volcano_plot))
  expect_identical(observed_volcano_plot$theme, expected_volcano_plot$theme)

  expect_silent(observed_volcano_plot <- plotVolcano(expected_contrasts_ann$BI_vs_AI, expected_parameters$results_directory, geneLabels=TRUE))
})

test_that("plotting volcano with bad parameters fails", {
  expect_error(plotVolcano(list()), "Could not plotVolcano: data is not a matrix, dataframe or DESeqResults object.")
  expect_error(plotVolcano(matrix()), "Could not plotVolcano: no results directory given.")
  
  bad_contrast <- expected_contrasts_ann$BI_vs_AI
  metadata(bad_contrast)$alpha <- 2
  expect_error(plotVolcano(bad_contrast, expected_parameters$results_directory), "Could not plotVolcano: alpha must be numeric, >0 and <1.")
  
  expect_error(plotVolcano(expected_contrasts_ann$BI_vs_AI, expected_parameters$results_directory, lfc='x'), "Could not plotVolcano: lfc is not numeric.")
  expect_error(plotVolcano(expected_contrasts_ann$BI_vs_AI, expected_parameters$results_directory, geneSymbols=c('x')), "Could not plotVolcano: length of geneSymbols is not the same as num rows in contrast.")
})

test_that("getting venn counts and plotting venn works", {
  expect_silent(observed_venn_counts <- getVennCounts(expected_contrasts[c(3,12)]))
  expect_identical(observed_venn_counts, expected_venn_counts)

  expect_silent(plotVenn(observed_venn_counts, expected_parameters$results_directory))
  venn_file <- file.path(expected_parameters$results_directory, 'images', 'DE_venn.png')
  expect_true(file.exists(venn_file))

  expect_silent(observed_venn_counts_ann <- getVennCounts(expected_contrasts_ann[c(3,12)]))
  expect_silent(plotVenn(observed_venn_counts_ann, expected_parameters$results_directory))
})

test_that("getting venn counts with >5 or <2 contrasts fails", {
  bad_contrast_list<- list(1,2,3,4,5,6)
  names(bad_contrast_list) <- c(1:6)
  expect_output(getVennCounts(bad_contrast_list), "A venn diagram could not be produced as there are more than 5 contrasts.")
  
  bad_contrast_list<- list(1)
  names(bad_contrast_list) <- 1
  expect_output(getVennCounts(bad_contrast_list), "A venn diagram could not be produced as there are not enough contrasts.")
})

test_that("getting and labelling top genes works", {
  gene_symbols <- paste0("test_", rownames(expected_contrasts[[3]]))
  expect_silent(observed_labelled_genes <- labelDEgenes(expected_contrasts[[3]], geneSymbol=gene_symbols, lfc=3, alpha=0.01))
  expect_silent(observed_top_genes <- getTopGenes(observed_labelled_genes, 10))
  expect_identical(observed_labelled_genes, expected_labelled_genes)
  expect_identical(observed_top_genes, expected_top_genes)
})
