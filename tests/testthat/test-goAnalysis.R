context("go analysis")

test_that("preparing GO table works", {
  expect_silent(observed_go_dt <- suppressMessages(prepareGOtable(expected_dds_de_ann, expected_contrasts_ann[["BI_vs_AI"]], 'BP')))
  expect_equal(observed_go_dt$x, expected_go_dt$x)
})

test_that("preparing GO data works", {
  expect_silent(observed_go_data <- suppressMessages(prepareGoData(expected_dds_de_ann, expected_contrasts_ann[["BI_vs_AI"]], 'BP')))
  expect_equal(observed_go_data, expected_go_data)
})

test_that("running GO analysis works", {
  expect_silent(observed_go_table <- suppressMessages(runGoAnalysis(expected_go_data)))
  expect_equal(observed_go_table, expected_go_table)
})

test_that("getting GO symbols works", {
  expect_silent(observed_go_symbols <- getGoSymbols(expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table))
  expect_equal(observed_go_symbols, expected_go_symbols)
})
