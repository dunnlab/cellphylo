#' subsample_combined_matrix
#'
#' @importFrom Seurat Read10X as.sparse
#' @importFrom dplyr filter
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to a combined matrix of cells from multiple species. Cell ids must be in cellphylo format.
#' @param sample_size The number of replicates per cell type group to be sampled
#' @param print Print out a file listing the cell ids of selected cells and the subsampled matrix in 10x format.
#'
#' @return mat.sparse The subsampled matrix in sparse matrix format.
#' @export
#'
#' @examples
subset_combined_matrix <- function(matrix_path, sample_size, print){

  #load matrix
  mat <- Read10X(matrix_path)
  #parse cell ids
  orig_ids <- colnames(mat)
  cluster_ids <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[3])})
  species_id <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[1])})
  species_cluster_ids <- paste0(species_id, "_", cluster_ids)

  #organize annotations into a dataframe
  ids.df <-data.frame(orig_ids = orig_ids, cluster_ids = cluster_ids, species_cluster_ids = species_cluster_ids)
  species_cluster_ids.unique <- unique(species_cluster_ids)
  length(species_cluster_ids.unique)


  subset_dfs <- lapply(species_cluster_ids.unique, function(x){

    ids.df.subset <- ids.df %>% filter(species_cluster_ids== x)
    random_subset <- sample(ids.df.subset$orig_ids, size = sample_size, replace=FALSE)
    random_subset <- as.data.frame(random_subset)
    return(random_subset)
  })

  combined_subset_dfs <- do.call(rbind, subset_dfs)
  n_subset <- length(species_cluster_ids.unique)*sample_size


  #Find the indices of the randomly selected cell barcodes in the orignal integrated matrix
  selected <- combined_subset_dfs$random_subset
  indices <- which(colnames(mat) %in% selected)

  if (print==TRUE){
    write.table(selected, paste0("aqhumor_cross_species_",n_subset, "_subset.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)}

  #subset and write out matrix
  mat.sub <-mat[,indices]
  mat.sparse <- as.sparse(mat.sub)

  if (print==TRUE){
    write10xCounts(paste0("aqhumor_cross_species_",n_subset, "_subset_mtx"), mat.sparse, version="3")
  }

  return(mat.sparse)

}
