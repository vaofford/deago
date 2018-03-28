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
  
  expect_silent(observed_go_dt_up <- suppressMessages(prepareGOtable(expected_go_tables[["BI_vs_AI_BP_up"]])))
  expect_equal(observed_go_dt_up$x$data, expected_go_dt_up$x$data, check.attributes = FALSE)
  
  expect_silent(observed_go_dt_down <- suppressMessages(prepareGOtable(expected_go_tables[["BI_vs_AI_BP_down"]])))
  expect_equal(observed_go_dt_down$x$data, expected_go_dt_down$x$data, check.attributes = FALSE)
})

test_that("preparing GO table gives message with no sig.genes", {
  expect_message(prepareGOtable(expected_go_tables[["AII_vs_AI_BP"]]), "A GO enrichment analysis could not be performed as there are no significantly differentially expressed genes for this contrast.")
})

test_that("getting gene identifiers works", {
  expect_silent(observed_go_identifiers <-getGOidentifiers (expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table))
  expect_equal(observed_go_identifiers, expected_go_identifiers)
  
  expect_silent(observed_go_identifiers_up <- getGOidentifiers (expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table, 'up'))
  expect_equal(observed_go_identifiers_up, expected_go_identifiers_up)
  
  expect_silent(observed_go_identifiers_down <- getGOidentifiers (expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table, 'down'))
  expect_equal(observed_go_identifiers_down, expected_go_identifiers_down)
})

test_that("getting gene symbols works", {
  expect_silent(observed_go_symbols <- getGOsymbols(expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table))
  expect_equal(observed_go_symbols, expected_go_symbols)
  
  expect_silent(observed_go_symbols_up <- getGOsymbols(expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table, 'up'))
  expect_equal(observed_go_symbols_up, expected_go_symbols_up)
  
  expect_silent(observed_go_symbols_down <- getGOsymbols(expected_contrasts_ann[["BI_vs_AI"]], expected_go_data, expected_go_table, 'down'))
  expect_equal(observed_go_symbols_down, expected_go_symbols_down)
})

test_that("writing GO tables works", {
  expect_silent(writeGOtables(expected_go_tables, expected_parameters$results_directory))
  expect_true(file.exists(file.path(expected_parameters$results_directory, "BI_vs_AI_BP.tsv")))
})



