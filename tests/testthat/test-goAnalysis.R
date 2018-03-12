context("go analysis")

test_that("running GO analysis works", {
  expect_silent(observed_go_tables <- suppressMessages(runGOanalysis(expected_dds_de_ann, list("BI_vs_AI"=expected_contrasts_ann[["BI_vs_AI"]]) , expected_parameters)))

#  observed_go_tables$BI_vs_AI_BP <- observed_go_tables$BI_vs_AI_BP[order(observed_go_tables$BI_vs_AI_BP$GO.ID, method='radix'),]
#  expected_go_tables$BI_vs_AI_BP <- expected_go_tables$BI_vs_AI_BP[order(expected_go_tables$BI_vs_AI_BP$GO.ID, method='radix'),]

#  expect_equal(observed_go_tables, expected_go_tables, check.attributes = FALSE)
})

test_that("preparing GO data works", {
  expect_silent(observed_go_data <- suppressMessages(prepareGOdata(expected_dds_de_ann, expected_contrasts_ann[["BI_vs_AI"]], c('BP'))))
  expect_equal(observed_go_data, expected_go_data, check.attributes = FALSE)
})

test_that("running topGO analysis works", {
  expect_silent(observed_go_table <- suppressMessages(topGOanalysis(expected_go_data)))

  # To satisfy R-devel differences in ordering
#  observed_go_table <- observed_go_table[order(observed_go_table$GO.ID, method='radix'),]
#  expected_go_table <- expected_go_table[order(expected_go_table$GO.ID, method='radix'),]

#  expect_equal(observed_go_table, expected_go_table, check.attributes = FALSE)
})

test_that("preparing GO table works", {
  expect_silent(observed_go_dt <- suppressMessages(prepareGOtable(expected_go_tables[["BI_vs_AI_BP"]])))
  expect_equal(observed_go_dt$x$data, expected_go_dt$x$data, check.attributes = FALSE)
})

test_that("getting GO symbols works", {
  expect_silent(observed_go_symbols <- getGOsymbols(expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table))
  expect_equal(observed_go_symbols, expected_go_symbols)
})

test_that("writing GO tables works", {
  expect_silent(writeGOtables(expected_go_tables, expected_parameters$results_directory))
  expect_true(file.exists(file.path(expected_parameters$results_directory, "BI_vs_AI_BP.tsv")))
})

