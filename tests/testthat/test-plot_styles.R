context("plot styles")

# Only general styles are tested here.  Plot-specific styles are tested in test_qcPlots.R or test_contrastPlots.R

test_that("theme_deago works", {
  tm <- theme_deago()
  expect_identical(tm$text$family, "serif")
  expect_identical(tm$text$size, 12)
  expect_identical(tm$axis.title$size, 14)
  expect_identical(tm$axis.title$face, 'bold')
  expect_identical(tm$axis.title$margin, margin(0, 15, 0, 0, 'pt'))
  expect_identical(tm$plot.title$size, 16)
  expect_identical(tm$plot.title$face, 'bold')
  expect_identical(tm$plot.title$hjust, 0.5)
})

test_that("theme_deago_legend works", {
  tm <- theme_deago_legend()
  expect_identical(tm$legend.title$size, 14)
  expect_identical(tm$legend.title$face, 'bold')
  expect_identical(tm$legend.title$vjust, 2)
  expect_identical(tm$legend.text$size, 12)
})


