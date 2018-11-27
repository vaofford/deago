test_directory <- tempdir()
test_images_directory <- file.path(test_directory, 'images')
if(!dir.exists(test_images_directory)) dir.create(test_images_directory)

stopifnot(dir.exists(test_directory))

load(file.path('deago-testdata.RData'))
expected_parameters$results_directory <- test_directory
expected_parameters$counts_directory <- file.path(expected_parameters$results_directory, expected_parameters$counts_directory)
expected_parameters$annotation_file <- file.path('deago-test-annotation.tsv')
if(!dir.exists(expected_parameters$counts_directory)) dir.create(expected_parameters$counts_directory)

bad_annotation_parameters <- list()
bad_annotation_parameters$annotation_file <- file.path('deago-test-bad-annotation.tsv')
switched_annotation_file <- file.path('deago-test-annotation-switch.tsv')

