context("config")

config_file <- file.path(expected_parameters$result_dir, 'test.config')

test_that("validating config works", {
  expect_silent(validateConfig(expected_parameters))

  expect_error(validateConfig(), "Could not validate config: need to specify list of parameters.")
  expect_error(validateConfig(c(1:3)), "Could not validate config: parameters are not a list.")

  essential_parameters = c('counts','target','result_dir')
  for (i in essential_parameters)
  {
    bad_parameters <- expected_parameters[ - which(names(expected_parameters) == i)]
    expect_error(validateConfig(bad_parameters), "Could not validate config: need to specify value for *")
  }

  bad_go_parameters <- c(expected_parameters[ - which(names(expected_parameters) == 'go_analysis' | names(expected_parameters) == 'annotation')], list('go_analysis'=1, 'annotation'=NULL))
  expect_error(validateConfig(bad_go_parameters), "Could not validate config: annotation not provided for go analysis")
})

test_that("writing config file works", {
  if(file.exists(config_file)) unlink(config_file)

  expect_error(writeConfig(), "Could not write config: need to specify file to write to.")
  expect_error(writeConfig(config_file), "Could not write config: need to specify list of parameters.")
  expect_error(writeConfig(config_file, parameters=c(1:3)), "Could not write config: parameters are not a list.")

  expect_silent(writeConfig(config_file, parameters=expected_parameters))
  expect_true(file.exists(config_file))
  expect_error(writeConfig(config_file, parameters=expected_parameters), "Could not write config: configuration file already exists.")
})

test_that("building config file works", {
  if(file.exists(config_file)) unlink(config_file)

  expect_error(buildConfig(), "Could not build config: need to specify file to write to.")
  expect_error(buildConfig(config_file), "Could not build config: need to specify list of parameters.")
  expect_error(buildConfig(config_file, parameters=c(1:3)), "Could not write config: parameters are not a list.")

  expect_silent(buildConfig(config_file, parameters=expected_parameters))
  expect_true(file.exists(config_file))
  expect_error(buildConfig(config_file, parameters=expected_parameters), "Could not build config: configuration file already exists.")
})

test_that("importing config file works", {
  if(!file.exists(config_file)) buildConfig(config_file, expected_parameters)

  expect_silent(observed_parameters <- importConfig(config_file))
  expect_identical(observed_parameters, expected_parameters)

  expect_error(importConfig(), "Could not import config: need to specify configuration file.")
  expect_error(importConfig('blah/blah/blah'), "Could not import config: configuration file does not exist.")

  empty_config_file <- file.path(expected_parameters$result_dir, 'empty_config')
  file.create(empty_config_file)
  expect_error(importConfig(empty_config_file), "Could not import config: configuration file is empty.")
  unlink(empty_config_file)

  bad_config_file <- file.path(expected_parameters$result_dir, 'bad_config')
  expect_silent(writeConfig(bad_config_file, parameters=list(1)))
  expect_error(importConfig(bad_config_file), "Could not import config: too few parameters.")
  unlink(bad_config_file)
})

