#' @title Read Configuration File
#'
#' @description \code{importConfig} reads configuration file into a list.
#'
#' @details \code{importConfig} reads a tab-delimited configuration file
#'   containing key/value pairs which represent the parameters for the analysis
#'   into a list.  The configuration file should contain two columns, the first
#'   containing the keys outlined below and the second, their values.
#'
#'   Required key/value pairs are: \describe{ \item{\strong{counts_directory}}{Directory
#'   where count data files are located.} \item{\strong{targets_file}}{Location of
#'   targets file containing the mappings between file names and experimental
#'   conditions.} \item{\strong{results_directory}}{Location where results directory
#'   is to be created.} }
#'
#'   Other accepted key/value pairs and their defaults are: \describe{
#'   \item{\strong{gene_ids}}{Name of column in count files which contains the
#'   gene identifiers [Default: geneID]} \item{\strong{alpha}}{Significance
#'   cut-off, see \link[DESeq2]{results} for more information [Default: 0.05]}
#'   \item{\strong{control}}{Name of condition in the targets file which was the
#'   control [Default: control]} \item{\strong{annotation_file}}{Location of
#'   annotation file, see \link[deago]{annotateDataset} [Default: NULL]}
#'   \item{\strong{keep_images}}{Give a value of 1 to keep images which are used
#'   in the HTML report, see \link[deago]{plotToFile} for more information
#'   [Default: 0] } \item{\strong{qc_only}}{Give a value of 0 to run DESeq2
#'   analysis functions or a value of 1 to report only the QC plots [Default:
#'   1]} \item{\strong{go_analysis}}{Give a value of 1 to run GO term enrichment
#'   analysis using \link{topGO}. If value of 1 given for go_analysis, an
#'   annotation file must be provided as a value for the annotation key
#'   [Default: 0]} }
#'
#'   Parameters are validated using \code{\link{validateConfig}}.
#'
#'   \code{importConfig} returns a list containing the parameter values with the
#'   names corresponding to their keys.
#'
#' @family configuration functions
#' @family import functions
#'
#' @param file A \link{character} string giving the name of the file.
#' @param path A \link{character} string giving the directory containing the
#'   file. This may be omitted if the file is in the current working directory.
#' @param sep A \link{character} used to separate the fields.
#'
#' @return A key/value \link{list} containing the parameters from the configuration file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  # Import tab-delimited config from working directory
#'  importConfig("config.txt")
#'
#'  # Import tab-delimited config from another directory
#'  importConfig("config.txt", path="/path/to/config")
#'
#'  # Import comma-delimited config from working directory
#'  importConfig("config.txt", sep=",")
#' }

importConfig <- function(file, path=NULL, sep="\t")
{
  if (missing(file)) stop("Could not import config: need to specify configuration file.")
  if (!is.null(path)) file <- file.path(path, file)

  if (!file.exists(file)) stop("Could not import config: configuration file does not exist.")
  if (file.info(file)$size == 0) stop("Could not import config: configuration file is empty.")

  parameters <- scan(file, what="list", sep="\n", quote = "\"", quiet=TRUE)
  if(length(parameters) < 2) stop("Could not import config: too few parameters.")

  parameters <- strsplit(parameters, sep)
  names(parameters) <- sapply(parameters, `[[`, 1)
  parameters <- lapply(parameters, `[`, -1)

  numeric_parameters <- c('alpha', 'keep_images', 'qc_only', 'go_analysis')
  for (i in numeric_parameters)
  {
    if(exists(i, where=parameters)) parameters[[i]] <- as.numeric(parameters[[i]])
  }

  parameters <- validateConfig(parameters)

  return(parameters)
}

#' @title Build Configuration File
#'
#' @description \code{buildConfig} writes a list of validated key/value pairs to
#'   a tab-delimited file.
#'
#' @details \code{buildConfig} takes a list of key/value pairs which represent
#'   the parameters for the analysis, validates them using
#'   \code{\link{validateConfig}} and then writes them to a tab-delimited
#'   configuration file using \code{\link{writeConfig}}.
#'
#'   See \code{\link{importConfig}} for required and accepted key/value pairs.
#'
#'   \code{buildConfig} returns the new configuration file location.
#'
#' @family configuration functions
#'
#' @param file A \link{character} string giving the name of the file.
#' @param path A \link{character} string giving the directory containing the
#'   file. This may be omitted if the file is in the current working directory.
#' @param parameters A \link{list} containing key/value pairs which define the
#'   parameters for the analysis.
#'
#' @return A \link{character} string containing the name of the output configuration file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  # Basic parameter list
#'  params <- list ( 'counts_directory' = "/path/to/counts.txt",
#'                   'targets_file'     = "/path/to/targets.txt",
#'                   'results_directory' = "/path/to/result_dir")
#'
#'  # Write tab-delimited config file to working directory
#'  buildConfig("config.txt", parameters=params)
#'
#'  # Write tab-delimited config file to another directory
#'  buildConfig("config.txt", path="/path/to/config", parameters=params)
#' }

buildConfig <- function(file, path=NULL, parameters)
{
  if (missing(file)) stop("Could not build config: need to specify file to write to.")
  if (!is.null(path)) file <- file.path(path, file)

  if (file.exists(file)) stop("Could not build config: configuration file already exists.")

  if (missing(parameters)) stop("Could not build config: need to specify list of parameters.")
  if (!is.list(parameters)) stop("Could not write config: parameters are not a list.")

  parameters <- validateConfig(parameters)
  writeConfig(file, path, parameters)

  return(file)
}

#' @title Write Configuration File
#' @description \code{writeConfig} writes collapsed list to a tab-delimited
#'   file.
#'
#' @details \code{writeConfig} collapses a list of key/value pairs representing
#'   the parameters for the analysis and writes them to a specified
#'   tab-delimited file.
#'
#'   \code{writeConfig} is used to write configuration files in
#'   \code{\link{buildConfig}}. See \code{\link{importConfig}} for required and
#'   accepted key/value pairs.
#'
#'   \code{writeConfig} returns the new configuration file location.
#'
#' @family configuration functions
#'
#' @param file character string: giving the name of the configuration file to
#'   write parameters.
#' @param path character string: giving the directory to write the file.
#' @param parameters list: giving key/value pairs of paramenters for analysis.
#'
#' @return A \link{character} string containing the name of the output file.
#'
#' @importFrom utils write.table
#' @export
#'
#' @examples
#' \dontrun{
#'  # Basic parameter list
#'  params <- list ( 'counts_directory' = "/path/to/counts.txt",
#'                   'targets_file'     = "/path/to/targets.txt",
#'                   'results_directory' = "/path/to/result_dir")
#'
#'  # Write tab-delimited config file to working directory
#'  writeConfig("config.txt", parameters=params)
#'
#'  # Write tab-delimited config file to another directory
#'  writeConfig("config.txt", path="/path/to/config", parameters=params)
#' }


writeConfig <- function(file, path=NULL, parameters)
{
  if (missing(file)) stop("Could not write config: need to specify file to write to.")
  if (!is.null(path)) file <- file.path(path, file)

  if (file.exists(file)) stop("Could not write config: configuration file already exists.")

  if (missing(parameters)) stop("Could not write config: need to specify list of parameters.")
  if (!is.list(parameters)) stop("Could not write config: parameters are not a list.")

  collapsed_parameters <- cbind(unlist(parameters))
  write.table(collapsed_parameters, file=file, sep="\t", col.names=FALSE)

  return(file)
}

#' @title Validate Parameter List
#'
#' @description \code{validateConfig} validates a list of parameter key/value
#'   pairs.
#'
#' @details \code{validateConfig} takes a list of key/value pairs which
#' represent the parameters for the analysis, validates them and then returns
#' the validated parameter key/value list.
#'
#' Required key/value pairs are: \describe{ \item{\strong{counts_directory}}{Directory
#' where count data files are located.} \item{\strong{targets_file}}{Location of
#' targets file which contains the mapping files and experimental conditions.}
#' \item{\strong{results_directory}}{Directory where results directory will be
#' created.} }
#'
#' Other accepted key/value pairs and their defaults are: \describe{
#' \item{\strong{gene_ids}}{Name of column in count files which contains the
#' gene identifiers [Default: geneID]} \item{\strong{alpha}}{[Significance
#' cut-off, see \link[DESeq2]{results} for more information Default: 0.05]}
#' \item{\strong{control}}{Name of condition in the targets file which was the
#' control [Default: control]} \item{\strong{annotation_file}}{Location of annotation
#' file, see \link[deago]{annotateDataset} [Default: NULL]}
#' \item{\strong{keep_images}}{Give a value of 1 to keep images which are used
#' in the HTML report, see \link[deago]{plotToFile} for more information
#' [Default: 0] } \item{\strong{qc_only}}{Give a value of 0 to run DESeq2
#' analysis functions or a value of 1 to report only the QC plots [Default: 1]}
#' \item{\strong{go_analysis}}{Give a value of 1 to run GO term enrichment
#' analysis using \link{topGO}. If value of 1 given for go_analysis, an
#' annotation file must be provided as a value for the annotation_file key [Default:
#' 0]} }
#'
#' \code{validateConfig} returns a list containing the validated parameter
#' values with the names corresponding to their keys.
#'
#' @family configuration functions
#' @family import functions
#'
#' @param parameters A \link{list} containing key/value pairs which define the
#'   parameters for the analysis (see \link[deago]{importConfig} for more
#'   information).
#'
#' @return A \link{list} containing the validated parameters.
#'
#' @importFrom utils write.table
#' @export
#'
#' @examples
#'  # Basic parameter list
#'  params <- list ( 'counts_directory' = "/path/to/counts.txt",
#'                   'targets_file'     = "/path/to/targets.txt",
#'                   'results_directory' = "/path/to/result_dir")
#'
#'  # Validate parameters
#'  validateConfig(params)


validateConfig <- function(parameters)
{
  if (missing(parameters)) stop("Could not validate config: need to specify list of parameters.")
  if (!is.list(parameters)) stop("Could not validate config: parameters are not a list.")

  essential_parameters = c('counts_directory','targets_file','results_directory')

  for (i in essential_parameters)
  {
    if(!exists(i, where=parameters)) stop(paste0("Could not validate config: need to specify value for ", i))
  }

  default_parameters <- list ( 'gene_ids'        = "geneID",
                               'alpha'           = 0.05,
                               'control'         = NULL,
                               'columns'         = 'condition',
                               'annotation_file' = NULL,
                               'keep_images'     = 0,
                               'qc_only'         = 1,
                               'go_analysis'     = 0)

  for (i in names(default_parameters))
  {
    if(!exists(i, where=parameters))
    {
      parameters[[i]] <- default_parameters[[i]]
    }
  }

  if (parameters$go_analysis == 1 && is.null(parameters$annotation_file)) stop("Could not validate config: annotation not provided for go analysis")
  if (parameters$qc_only == 1) parameters$go_analysis <- 0

  return(parameters)
}

