#' wrangle_matrix
#'
#' Converts cell identifiers from van Zyl et al. (2020)'s format to cellphylo format.
#'
#' @importFrom Seurat Read10X
#' @importFrom DropletUtils write10xCounts
#' @importFrom utils read.table
#' @importFrom dplyr mutate left_join
#'
#' @param path_to_mat Path to matrix to be wrangled into cellphylo format. Matrix must be in 10x format.
#' @param metadata_file Path to wrangled metadata file from van Zyl et al. (2020)
#' @param species Name of species the matrix belongs to
#'
#' @return mat A wrangled matrix with cell ids in cellphylo format
#' @export
#'
#' @examples
wrangle_matrix <- function(path_to_mat, metadata_file, species){

  #load wrangled meta data file
  meta_data <- read.table(file=metadata_file, header=TRUE, sep="\t")

  #load matrix
  mat <- Read10X(path_to_mat)
  #extract cell ids in order they are found in matrix
  matrix_cell_ids <- colnames(mat) %>% as.data.frame()
  colnames(matrix_cell_ids) <- "cell_identifier"

  #reorder meta data df to match mat colnames
  names.df <- left_join(matrix_cell_ids, meta_data, by= "cell_identifier")

  #parse out meta data info
  sample_names <- names.df$sample_id
  cell_barcode_ids <- names.df$cell_barcode
  cluster_id <- names.df$cluster_id
  cell_type_id <- gsub('[0-9]+', "",names.df$cluster_id)
  species_id <- species

  #create new cell names
  names.df <- names.df %>% mutate(new_cell_ids = paste0(species_id, "_", sample_names,"_",cluster_id,"_", cell_type_id,"_", cell_barcode_ids))

  #replace matrix names
  colnames(mat) <- names.df$new_cell_ids

  #write the new matrix
  write10xCounts(paste0("aqhumor_UMI_wrangled_",species), mat, version="3")

  if(!dir.exists("matrix/wrangled")){
    dir.create("matrix/wrangled")
  }

  file.rename(paste0("aqhumor_UMI_wrangled_",species), paste0("matrix/wrangled/aqhumor_UMI_wrangled_",species))

  return(mat)

}
