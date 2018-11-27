context("results directory")

test_that("making results and image directories works", {
  results_dir <- makeResultDir(path=expected_parameters$results_directory, keep_images = 1)
  expect_true(dir.exists(results_dir))
  expect_silent(checkDirExists(results_dir))

  image_dir <- file.path(results_dir,"images")
  expect_true(dir.exists(image_dir))
  expect_silent(checkDirExists(image_dir))

  expect_error(makeResultDir(), "Could not create results directory: must specify output directory.")
  expect_error(makeResultDir(path='blah/blah/blah'), "Could not create results directory: specified path does not exist")
  expect_error(checkDirExists('dirDoesNotExist'), "dirDoesNotExist does not exist")
})
