#' convert_to_10x
#'
#' @importFrom utils read.csv
#' @importFrom R.utils gunzip
#' @importFrom tibble column_to_rownames
#' @importFrom DropletUtils write10xCounts
#' @importFrom Seurat as.sparse
#'
#' @param path_to_csv Path to csv file
#' @param dir_name Name for 10x matrix directory
#'
#' @return NULL a directory with a matrix in 10x format (version 3)
#' @export
#'
#' @examples
#'
convert_to_10x <- function(path_to_csv, dir_name){

  #gunzip
  unzipped <- gunzip(path_to_csv, remove=FALSE)
  #read csv
  csv <- read.csv(unzipped)
  #make rownames gene ids
  csv <- column_to_rownames(csv, var = "X")
  #convert to sparse
  sparse <- as.sparse(csv)
  #write new 10x matrix
  write10xCounts(paste0(dir_name), sparse, version="3")
  #remove expanded file
  remove_file <- gsub(".gz", "", path_to_csv)
  file.remove(remove_file)
} #close fun

