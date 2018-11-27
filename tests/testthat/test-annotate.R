context("annotation")

test_that("importing annotation and annotating DESeqDataSet works", {
  expect_error(annotateDataset(), 'argument "dds" is missing, with no default')
  expect_error(annotateDataset(expected_dds_de), 'argument "parameters" is missing, with no default')
  expect_error(annotateDataset(expected_dds_de, list()), "Could not import annotation: no annotation file provided.")
  expect_error(annotateDataset(expected_dds_de, bad_annotation_parameters), "Could not import annotation: no dataset identifiers found.")
  expect_silent(observed_dds_de_ann <- annotateDataset(expected_dds_de, expected_parameters))
  expect_identical(observed_dds_de_ann, expected_dds_de_ann)
})

test_that("importing annotation works", {
  expect_error(importAnnotation(), 'argument "file" is missing, with no default')
  expect_error(importAnnotation('noFile'), 'noFile does not exist')
  expect_silent(observed_annotation <- importAnnotation(expected_parameters$annotation_file))
  expect_identical(observed_annotation, expected_annotation)
  expect_silent(observed_annotation <- importAnnotation(switched_annotation_file))
})

test_that("importing annotation with too few columns fails", {
  bad_ann_file <- file.path(expected_parameters$results_directory, "too_few.txt")
  write.table(as.data.frame(matrix(0, ncol = 1, nrow = 5)), bad_ann_file, sep="\t", row.names = FALSE, col.names=FALSE)
  expect_error(importAnnotation(bad_ann_file), "Could not import annotation: too few columns in annotation file.")
  unlink(bad_ann_file)
})

test_that("importing annotation with too many columns fails", {
  bad_ann_file <- file.path(expected_parameters$results_directory, "too_few.txt")
  write.table(as.data.frame(matrix(0, ncol = 4, nrow = 5)), bad_ann_file, sep="\t", row.names = FALSE, col.names=FALSE)
  expect_error(importAnnotation(bad_ann_file), "Could not import annotation: too many columns in annotation file.")
  unlink(bad_ann_file)
})

test_that("adding annotation with GO terms in multiple columns fails", {
  bad_annotation <- cbind(expected_annotation[,c(1,3)], expected_annotation[,3])
  expect_error(addAnnotationsToDataSet(expected_dds_de, bad_annotation), "Could not import annotation: GO terms were found in multiple columns.")
})

test_that("adding annotation with 3 columns and no GO terms fails", {
  bad_annotation <- cbind(expected_annotation[,1:2], c(1:1000))
  expect_error(addAnnotationsToDataSet(expected_dds_de, bad_annotation), "Could not import annotation: no GO terms identified.")
})

test_that("adding annotation with wrong column number fails", {
  bad_annotation <- data.frame('x'=1:5)
  expect_error(addAnnotationsToDataSet(expected_dds_de, bad_annotation), "Could not import annotation: too few columns in annotation.")
  
  bad_annotation <- data.frame('a'=1:5,'b'=1:5,'c'=1:5,'d'=1:5)
  expect_error(addAnnotationsToDataSet(expected_dds_de, bad_annotation), "Could not import annotation: too many columns in annotation.")
})

test_that("adding annotation with genes only works", {
  expected_dds_de_ann_geneOnly <- expected_dds_de_ann
  metadata(expected_dds_de_ann_geneOnly) <- metadata(expected_dds_de_ann_geneOnly)[1]

  expect_silent(observed_dds_de_ann <- addAnnotationsToDataSet(expected_dds_de, expected_annotation[,1:2]))
  expect_equal(observed_dds_de_ann, expected_dds_de_ann_geneOnly, check.attributes = FALSE)
})

test_that("adding annotation with genes only works", {
  expected_dds_de_ann_goOnly <- expected_dds_de_ann
  mcols(expected_dds_de_ann_goOnly) <- mcols(expected_dds_de_ann_goOnly)[1:40]

  expect_silent(observed_dds_de_ann <- addAnnotationsToDataSet(expected_dds_de, expected_annotation[,c(1,3)]))
  expect_equal(observed_dds_de_ann, expected_dds_de_ann_goOnly, check.attributes = FALSE)
})

test_that("parsing go terms works", {
  expect_silent(observed_go_list <- getGOlist(expected_annotation, 3))
  expect_identical(observed_go_list, expected_go_list)

  expect_silent(observed_go_list <- getGOlist(expected_annotation[,c(1,3)], 2))
  expect_identical(observed_go_list, expected_go_list)

  expect_error(getGOlist(expected_annotation, 2), "Could not import annotation: no GO terms identified.")
  expect_error(getGOlist(expected_annotation[1,], 3), "Could not import annotation: GO list not created.")
})

test_that("parsing gene names works", {
  expect_silent(observed_gene_list <- getGeneSymbols(expected_dds_de, expected_annotation, 2)) # go terms and gene names
  expect_identical(observed_gene_list, expected_gene_list)

  expect_silent(observed_gene_list <- getGeneSymbols(expected_dds_de, expected_annotation[,1:2], 2))
  expect_identical(observed_gene_list, expected_gene_list)

  expect_error(getGeneSymbols(expected_dds_de, expected_annotation[,2:3], 2), "Could not import annotation: no dataset identifiers found.")
})
