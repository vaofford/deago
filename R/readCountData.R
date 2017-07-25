#' @title Read Count Data
#'
#' @description Generates a count matrix from a list of files containing count
#'   data.
#'
#' @details \code{readCountData} is the main data import function for
#'   \code{deago}.  Before importing count data, target information should be
#'   imported using \code{\link{importTargets}} which will map the count data
#'   files to the experimental conditions applied for each sample.  The targets
#'   dataframe should contain the following columns:
#'
#'   \describe{ \item{\strong{filename}}{name of the file containing the count
#'   data for the sample} \item{\strong{condition}}{experimental treatment or
#'   condition that was applied} \item{\strong{replicate}}{replicate identifier
#'   - can be numeric or character} \item{\strong{label}}{unique sample
#'   identifier comprised of the \code{condition} and \code{replicate}}}
#'
#'   For more information on the targets file and dataframe see
#'   \code{\link{importTargets}}.
#'
#'   When importing from individual files, count data will be imported using the
#'   \code{filename} and \code{label} columns in the targets dataframe. Each row
#'   in the targets file represents a sample, each sample in the targets file
#'   has a file name and a unique label.  The column containing the gene
#'   identifiers must have the same column name in each of the count data files.
#'   This column name should be specified using \code{id_column}. The order of
#'   these gene identifiers should be the same across all of the count data
#'   files being imported.
#'
#'   By default, \code{deago} assumes that the count data column is in the final
#'   column in the file. If this is not the case, a column number must be
#'   specified using \code{data_column}.  The reason that \code{deago} uses a
#'   column number instead of a column name is because many of the different
#'   count data programs use a filename as the header for the count column which
#'   would differ between files. As a file may have been renamed between
#'   creation and analysis, it is also simpler to make sure that the counts are
#'   in the same place in each of the files.
#'
#'   Once imported, the first column of the \code{deago} count matrix returned
#'   by \code{readCountData} will contain the gene identifiers. The remaining
#'   columns will contain the count data for each sample with the column name
#'   being the unique label associated with that sample in the targets
#'   dataframe.
#'
#' @family count data functions
#' @family import functions
#'
#' @param targets  A \link{data.frame} containing the target information mapping
#'   samples to their count data files and experimental conditions.  See
#'   \code{\link{importTargets}} for more information.
#' @param path A \link{character} string giving the directory containing the
#'   file. This may be omitted if the file is in the current working directory.
#' @param id_column A \link{character} string giving the name of the column
#'   which contains the gene identifiers, should be identical for all of the
#'   count data files.
#' @param data_column A \link{numeric} which is the index of the column in the
#'   files which contains the count data
#' @param skip A \link{numeric} which is the number of non-column-header lines
#'   to skip at start of the file
#' @param sep A \link{character} which is field separator
#'
#' @return dataframe: containing count data
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' readCountData(targets, id_column="genes", data_column=2)
#' }

readCountData <- function (targets, path, id_column, data_column, skip=0, sep="\t")
{
  if (missing(targets)) stop("Could not import counts: need to specify targets dataframe.")
  if (!is.numeric(data_column)) stop("Could not import counts: data column is not numeric.")
  if (!is.numeric(skip)) stop("Could not import counts: skip is not numeric.")

  if (!"filename" %in% colnames(targets)) stop("Could not import counts: filename column not found.")
  if (!"label" %in% colnames(targets)) stop("Could not import counts: label column not found.")

  file_count <- 1
  all_counts <- list()

  for (file in targets$filename)
  {
    if(!is.null(path)) file <- file.path(path, file)
    imported_counts <- read.table(file, header=TRUE, skip=skip, check.names=FALSE, stringsAsFactors = FALSE, sep=sep)

    if (!id_column %in% colnames(imported_counts)) stop(paste0("Could not import counts: gene id column ", id_column, "not found in ", file))

    current_counts <- suppressWarnings(imported_counts[order(imported_counts[,id_column]),])

    if (file_count == 1) gene_ids <- current_counts[,id_column]
    if (!identical(gene_ids, current_counts[,id_column])) stop(paste0("Could not import counts: gene id order differs in ", file))

    if (!missing(data_column) && data_column > ncol(current_counts)) stop(paste0("Could not import counts: data column greater than column number in ", file))
    if (missing(data_column)) data_column <- ncol(current_counts)

    all_counts[[file_count]] <- current_counts[,data_column]

    file_count <- file_count + 1
  }

  count_matrix <- do.call(cbind,all_counts)
  count_matrix <- as.data.frame(count_matrix, row.names=gene_ids)
  names(count_matrix) <- targets$label

  if (nrow(count_matrix) != length(gene_ids)) stop (paste0("Could not import counts: number of rows in count list( ",nrow(count_matrix)," ) not equal to number of genes in input files ",length(gene_ids)))
  if (ncol(count_matrix) != length(targets$filename)) stop (paste0("Could not import counts: length of count list( ",ncol(count_matrix)," ) not equal to number of filenames provided ",length(targets$filename)))

  count_matrix <- as.matrix(count_matrix)
  return(count_matrix)
}



