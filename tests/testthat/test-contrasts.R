context("contrasts")

test_that("getting contrasts works", {
  expect_silent(observed_contrasts <- getContrasts(expected_dds_de, expected_parameters))
  expect_is(observed_contrasts[['BI_vs_AI']], 'DESeqResults')
  expect_equal(observed_contrasts, expected_contrasts, check.attributes = FALSE)

  expect_silent(observed_contrasts_ann <- getContrasts(expected_dds_de_ann, expected_parameters))
  expect_is(observed_contrasts_ann[['BI_vs_AI']], 'DESeqResults')
  expect_equal(observed_contrasts_ann, expected_contrasts_ann, check.attributes = FALSE)
})

test_that("preparing contrast works", {
  expect_silent(observed_prepared_contrast <- prepareContrast(expected_dds_de_ann, expected_contrasts_ann[["BI_vs_AI"]]))
  expect_identical(observed_prepared_contrast, expected_prepared_contrast)
})

test_that("writing contrast works", {
  expect_silent(writeContrasts(expected_dds_de_ann, expected_contrasts_ann, expected_parameters$results_directory))
  for (i in names(expected_contrasts))
  {
    expect_true(file.exists(file.path(expected_parameters$results_directory, paste0(i,"_q0.05.txt"))))
  }
})

test_that("summarising contrasts works", {
  expect_silent(observed_contrast_summary <- contrastSummary(expected_contrasts_ann, expected_parameters))
  expect_is(observed_contrast_summary, 'data.frame')
  expect_identical(observed_contrast_summary, expected_contrast_summary)
})

test_that("preparing contrast table works", {
  expect_silent(observed_contrast_table <- prepareContrastTable(expected_contrasts_ann[["BI_vs_AI"]]))

  expect_is(observed_contrast_table, 'datatables')
  expect_is(observed_contrast_table, 'htmlwidget')

  # To satisfy R-devel differences in ordering
  observed_contrast_table$x$data <- observed_contrast_table$x$data[order(row.names(observed_contrast_table$x$data), method='radix'),]
  expected_contrast_table$x$data <- expected_contrast_table$x$data[order(row.names(expected_contrast_table$x$data), method='radix'),]

  expect_equal(observed_contrast_table$x$data, expected_contrast_table$x$data, check.attributes = FALSE)
})

test_that("preparing empty contrast table gives message", {
  expect_message(prepareContrastTable(expected_contrasts_ann[["AII_vs_AI"]]), "Could not generate contrast table for*")
})
