context("targets")

test_that("importing targets file works", {
  expect_error(importTargets(), "Could not import targets: need to specify targets file.")

  targets_file <- file.path(expected_parameters$results_directory, "test_targets.txt")
  write.table(target_df, file=targets_file, quote=FALSE, row.names=FALSE, sep="\t")
  expect_true(file.exists(targets_file))

  expect_silent(observed_targets <- importTargets(expected_parameters$target, path=expected_parameters$results_directory))
  expect_is(observed_targets, 'data.frame')
  expect_identical(observed_targets, expected_targets)
})

test_that ("validating targets file works", {
  expect_true(validateTargets(expected_targets[,1:3]))

  expect_error(validateTargets(),"Could not import targets: need to specify targets dataframe.")
  expect_error(validateTargets(expected_targets[,1:2]), "Could not import targets: too few columns in targets dataframe.")

  expected_columns <- c('filename', 'replicate', 'condition')
  for (coln in expected_columns)
  {
    bad_dummy_targets <- subset(expected_targets, select = names(expected_targets) != coln)
    expect_error(validateTargets(bad_dummy_targets), "* column does not exist")
  }
})

