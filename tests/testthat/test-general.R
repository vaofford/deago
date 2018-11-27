context("general functions")

test_that("rounding to base works", {
  expect_equal(roundToBase(17, 2, 'up'), 18)
  expect_equal(roundToBase(17, 10, 'up'), 20)
  expect_equal(roundToBase(17, 10, 'down'), 10)
  expect_equal(roundToBase(-17, 10, 'down'), -20)
})

test_that("rounding to base fails", {
  expect_error(roundToBase(17, 2, 'bad'),"Could not roundToBase, direction not up or down.")
})

test_that("rounding to factor of ten works", {
  expect_equal(roundToFactorOfTen(19), 20)
  expect_equal(roundToFactorOfTen(194), 200)
  expect_equal(roundToFactorOfTen(1094), 2000)
})

test_that("rounding to factor of ten fails", {
  expect_error(roundToFactorOfTen(c(19,4)), "'x' must be of length 1")
})
