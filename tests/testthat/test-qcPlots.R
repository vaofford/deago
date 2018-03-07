context("qc plots")

test_that("plotting read counts works", {
  testPlot('plotReadCounts', expected_dds_de, expected_rc_plot, 'reads_per_sample.png', expected_parameters$results_directory)
})

test_that("plotting null counts works", {
  testPlot('plotNullCounts', expected_dds_de, expected_nc_plot, 'null_counts_per_sample.png', expected_parameters$results_directory)
})

test_that("plotting sample distances works", {
  testPlot('plotSampleDistances', expected_dds_de, expected_sd_plot, 'sample_distances.png', expected_parameters$results_directory)
})

test_that("preparing PC list works", {
  expect_silent(observed_pc_list <- getPrincipalComponents(expected_dds_de))
  expect_equal(observed_pc_list, expected_pc_list)
})

test_that("plotting PCA works", {
  testPlot('pcaPlot', expected_pc_list, expected_pca_plot, 'pca.png', expected_parameters$results_directory)
})

test_that("plotting PCA scree works", {
  testPlot('pcaScreePlot', expected_pc_list, expected_pca_scree_plot, 'pca_scree.png', expected_parameters$results_directory)
})

test_that("preparing PCA table works", {
  expect_silent(observed_pca_table <- pcaSummary(expected_pc_list))
  expect_is(observed_pca_table, 'datatables')
  expect_is(observed_pca_table, 'htmlwidget')
  expect_equal(observed_pca_table$data, expected_pca_table$data)
})

test_that("plotting Cook's distances works", {
  testPlot('plotCooks', expected_dds_de, expected_cooks_plot, 'cooks_distances.png', expected_parameters$results_directory, suppressPlotMessages=TRUE)
})

test_that("plotting normalised densities works", {
  testPlot('plotDensity', expected_dds_de, expected_density_plot, 'normalised_density.png', expected_parameters$results_directory)
})

test_that("plotting dispersion estimates works", {
  testPlot('plotDispersionEstimates', expected_dds_de, expected_dispersion_plot, 'dispersion_estimates.png', expected_parameters$results_directory)
})

