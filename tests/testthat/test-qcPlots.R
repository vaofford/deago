context("qc plots")

test_that("plotting read counts works", {
  testPlot('plotReadCounts', expected_dds_de, expected_rc_plot, 'reads_per_sample.png', expected_parameters$results_directory)
})

test_that("plotting sample distances works", {
  testPlot('plotSampleDistances', expected_dds_de, expected_sd_plot, 'sample_distances.png', expected_parameters$results_directory)
})

test_that("plotting PCA works", {
  testPlot('pcaPlot', expected_dds_de, expected_pca_plot, 'pca.png', expected_parameters$results_directory)
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

