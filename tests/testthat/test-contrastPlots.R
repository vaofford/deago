context("contrast plots")

test_that("plotting MA works", {
  expect_silent(observed_ma_plot <- plotContrastMA(expected_contrasts$BI_vs_AI, expected_parameters$result_dir, geneLabels=TRUE))
  expect_is(observed_ma_plot, 'ggplot')

  observed_ma_file <- file.path(expected_parameters$result_dir, 'images', 'BI_vs_AI_MA.png')
  expect_true(file.exists(observed_ma_file))

  expect_identical(layer_data(observed_ma_plot), layer_data(expected_ma_plot))
  expect_equal(layer_scales(observed_ma_plot), layer_scales(expected_ma_plot))
  expect_identical(observed_ma_plot$theme, expected_ma_plot$theme)

  expect_silent(observed_ma_plot <- plotContrastMA(expected_contrasts_ann$BI_vs_AI, expected_parameters$result_dir, geneLabels=TRUE))
})

test_that("plotting volcano works", {
  expect_silent(observed_volcano_plot <- plotVolcano(expected_contrasts$BI_vs_AI, expected_parameters$result_dir, geneLabels=TRUE))
  expect_is(observed_volcano_plot, 'ggplot')

  observed_volcano_file <- file.path(expected_parameters$result_dir, 'images', 'BI_vs_AI_volcano.png')
  expect_true(file.exists(observed_volcano_file))

  expect_equal(layer_data(observed_volcano_plot), layer_data(expected_volcano_plot))
  expect_equal(layer_scales(observed_volcano_plot), layer_scales(expected_volcano_plot))
  expect_identical(observed_volcano_plot$theme, expected_volcano_plot$theme)

  expect_silent(observed_volcano_plot <- plotVolcano(expected_contrasts_ann$BI_vs_AI, expected_parameters$result_dir, geneLabels=TRUE))
})

test_that("getting venn counts and plotting venn works", {
  expect_silent(observed_venn_counts <- getVennCounts(expected_contrasts[c(3,12)]))
  expect_identical(observed_venn_counts, expected_venn_counts)

  expect_silent(plotVenn(observed_venn_counts, expected_parameters$result_dir))
  venn_file <- file.path(expected_parameters$result_dir, 'images', 'DE_venn.png')
  expect_true(file.exists(venn_file))

  expect_silent(observed_venn_counts_ann <- getVennCounts(expected_contrasts_ann[c(3,12)]))
  expect_silent(plotVenn(observed_venn_counts_ann, expected_parameters$result_dir))
})

test_that("getting and labelling top genes works", {
  gene_symbols <- paste0("test_", rownames(expected_contrasts[[3]]))
  expect_silent(observed_labelled_genes <- labelDEgenes(expected_contrasts[[3]], geneSymbol=gene_symbols, lfc=3, alpha=0.01))
  expect_silent(observed_top_genes <- getTopGenes(observed_labelled_genes, 10))
  expect_identical(observed_labelled_genes, expected_labelled_genes)
  expect_identical(observed_top_genes, expected_top_genes)
})
