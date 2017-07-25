test_directory <- tempdir()
test_images_directory <- file.path(test_directory, 'images')
if(!dir.exists(test_images_directory)) dir.create(test_images_directory)

stopifnot(dir.exists(test_directory))

load(file.path('deago-testdata.RData'))
expected_parameters$result_dir <- test_directory
expected_parameters$counts <-file.path(expected_parameters$result_dir, expected_parameters$counts)
expected_parameters$annotation <-file.path('deago-test-annotation.tsv')
if(!dir.exists(expected_parameters$counts)) dir.create(expected_parameters$counts)



