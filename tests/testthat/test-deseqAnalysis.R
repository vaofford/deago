context("DESeq analysis")

expected_coldata <- as.data.frame(colData(expected_dds))
expected_countdata <- counts(expected_dds)

test_that("preparing coldata works", {
  expect_silent(observed_coldata <- prepareColData(expected_targets, expected_countdata, expected_parameters))
  expect_equal(observed_coldata, expected_coldata)

  relevelled_coldata <- expected_coldata
  relevelled_coldata$condition <- relevel(relevelled_coldata$condition, ref='bi')
  expect_silent(observed_levelled_coldata <- prepareColData(expected_targets, expected_countdata, list('control'='bi')))
  expect_identical(observed_levelled_coldata, relevelled_coldata)
})

test_that("constructing DESeqDataSet works", {
  expect_silent(observed_dds <- constructDESeqDataset(expected_countdata, expected_coldata))
  expect_is(observed_dds, 'DESeqDataSet')
  expect_equivalent(observed_dds, expected_dds)
})

test_that("running DESeq analysis works", {
  expect_silent(observed_dds_de <- suppressMessages(runDESeqAnalysis(expected_targets, expected_countdata, expected_parameters)))
  expect_is(observed_dds_de, 'DESeqDataSet')
  expect_equivalent(observed_dds_de, expected_dds_de)
})
