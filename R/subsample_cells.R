#' subsample_cells
#' @importFrom Seurat Read10X as.sparse
#' @importFrom dplyr filter
#' @importFrom utils write.table
#' @importFrom DropletUtils write10xCounts
#'
#' @param path_to_matrix Path to matrix you want to subsample from. Cell ids must be in cellphylo format.
#' @param n_reps Number of replicates per each cell type group to be in the subsampled matrix.
#'
#' @return mat.sub A matrix created by random subsampling of the input matrix. All cell ids must be in cellphylo format.
#'
#' @examples
subsample_cells <- function(path_to_matrix, n_reps){

  print("Reading in matrix...")
  mat <- Read10X(path_to_matrix)
  print("Matrix loaded.")

  #species identity
  species <- sapply(strsplit(colnames(mat)[1], "_"), function(v){return(v[1])})
  #original full cell id
  orig_ids <- colnames(mat)
  #cluster ids (cluster num + cell type)
  cluster_ids <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[3])})
  #organize ids into a dataframe
  ids.df <- data.frame(orig_ids = orig_ids, cluster_id = cluster_ids)

  #how many cluster ids
  cluster_ids.unique <- unique(cluster_ids)

  subset_dfs <- lapply(cluster_ids.unique, function(x){

    ids.df.subset <- ids.df %>% filter(cluster_id == x)
    random_subset <- sample(ids.df.subset$orig_ids, size = n_reps, replace=FALSE)
    random_subset <- as.data.frame(random_subset)
    return(random_subset)
  })

  combined_subset_dfs <- do.call(rbind, subset_dfs)


  #Find the cell ids of the randomly selected cell barcodes in the orignal integrated matrix. Old approach used indices. Save ids.
  selected <- combined_subset_dfs$random_subset
  write.table(selected, paste0("subset_cells_integrated_", species, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

  #subset matrix
  indices <- which(colnames(mat) %in% selected)
  mat.sub <-mat[,indices]

  #write matrix
  mat.sparse <- as.sparse(mat.sub)
  write10xCounts(paste0("integrated_subset_mtx_", species), mat.sub, version="3")


  return(mat.sub)

} #close function
