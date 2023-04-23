#' convert_to_10x
#'
#' Converts a UMI count matrix in CSV format to 10x-style feature-barcode matrices.
#'
#' @importFrom utils read.csv
#' @importFrom R.utils gunzip
#' @importFrom tibble column_to_rownames
#' @importFrom DropletUtils write10xCounts
#' @importFrom Seurat as.sparse
#'
#' @param path_to_csv Path to csv file
#' @param dir_name Name for the new 10x matrix directory
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
  write10xCounts(dir_name, sparse, version="3")
  #create the 10x folder if not already present
  if(!dir.exists("matrix/10x")){
    dir.create("matrix/10x")
  }
  #move this file into the 10x folder
  file.rename(dir_name, paste0("matrix/10x/",dir_name))
  #remove expanded file
  remove_file <- gsub(".gz", "", path_to_csv)
  file.remove(remove_file)
} #close fun

