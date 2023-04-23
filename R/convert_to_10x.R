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
#' @param species_name Name of species matrix belongs to.
#'
#' @return sparse The matrix in sparse matrix format and a directory with a matrix in 10x feature-barcode format (version 3)
#' @export
#'
#' @examples
#'
convert_to_10x <- function(path_to_csv, species_name){

  #gunzip
  unzipped <- gunzip(path_to_csv, remove=FALSE)
  #read csv
  csv <- read.csv(unzipped)
  #make rownames gene ids
  csv <- column_to_rownames(csv, var = "X")
  #convert to sparse
  sparse <- as.sparse(csv)
  #write new 10x matrix
  write10xCounts(paste0("matrix_10x_", species_name), sparse, version="3")
  #create a folder for matrices, if not already present
  if(!dir.exists("matrix")){
    dir.create("matrix")
  }
  #create the 10x folder if not already present
  if(!dir.exists("matrix/10x")){
    dir.create("matrix/10x")
  }
  #move this file into the 10x folder
  file.rename(paste0("matrix_10x_", species_name), paste0("matrix/10x/matrix_10x_", species_name))
  #remove expanded file
  remove_file <- gsub(".gz", "", path_to_csv)
  file.remove(remove_file)

  return(sparse)
} #close fun

