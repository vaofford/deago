context("count data")

count_dir <- expected_parameters$counts_directory

for (i in colnames(expected_dds))
{
  count_file <- file.path(count_dir, paste0(i, '.tsv'))
  count_data <- data.frame('gene_names'=rownames(expected_dds), 'counts'=as.vector(counts(expected_dds)[,i]))
  write.table(count_data, file=count_file, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
}

test_that("reading in count data works", {
  expected_counts <- counts(expected_dds)
  expected_counts <- expected_counts[order(row.names(expected_counts)),]
  colnames(expected_counts) <- expected_targets$label
  observed_counts <- readCountData(expected_targets, path=count_dir, id_column='gene_names', data_column=2)
  expect_is(observed_counts, 'matrix')
  expect_equal(expected_counts, observed_counts)
})

test_that("missing parameters fails", {
  expect_error(readCountData(expected_targets, path=count_dir), "* is missing, with no default")
  expect_error(readCountData(path=count_dir, id_column='gene_names', data_column=2), "Could not import counts: need to specify targets dataframe.")
  expect_error(readCountData(expected_targets, path=count_dir, id_column='gene_names', data_column='x'), "Could not import counts: data column is not numeric.")
  expect_error(readCountData(expected_targets, path=count_dir, id_column='gene_names', data_column=2, skip='x'), "Could not import counts: skip is not numeric.")
})

test_that("bad targets fails", {
  expect_error(readCountData(expected_targets[,1:3], path=count_dir, id_column='gene_names', data_column=2), "Could not import counts: label column not found.")
  expect_error(readCountData(expected_targets[,1:2], path=count_dir, id_column='gene_names', data_column=2), "Could not import counts: filename column not found.")
})

test_that("bad parameters fails", {
  expect_error(readCountData(expected_targets, path=count_dir, id_column='gene_name', data_column=2), "Could not import counts: gene id column *")
  expect_error(readCountData(expected_targets, path=count_dir, id_column='gene_names', data_column=3), "Could not import counts: data column greater than column number in *")
})

test_that("different gene id order fails", {
  bad_targets <- rbind(expected_targets, data.frame('condition'='AI', 'replicate'='4', 'filename'='bad_sample.tsv', 'label'='ai_4'))
  bad_count_data <- data.frame('gene_names'=rownames(expected_dds), 'counts'=as.vector(counts(expected_dds)[,1]))
  bad_count_data <- bad_count_data[c(1,3,2),]
  bad_count_file <- file.path(count_dir, 'bad_sample.tsv')
  write.table(bad_count_data, file=bad_count_file, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
})

unlink(count_dir)
