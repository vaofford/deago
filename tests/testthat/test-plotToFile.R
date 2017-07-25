context("plotting to file")

test_that("plotting to file works", {
  test_plot <- ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point()
  expect_silent(p <- plotToFile(test_plot, expected_parameters$result_dir, 'test_plot.png', width=800, height=800))
  expect_true(file.exists(file.path(expected_parameters$result_dir, 'images', 'test_plot.png')))

  expect_error(plotToFile(), "Could not write plot: no plot given.")
  expect_error(plotToFile(test_plot), "Could not write plot: no results directory given.")
  expect_error(plotToFile(test_plot, expected_parameters$result_dir), "Could not write plot: no filename given.")
})
