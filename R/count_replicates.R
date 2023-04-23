#' count_replicates
#'
#' This function counts the number of replicates per cell type group.
#'
#' @importFrom Seurat Read10X
#'
#'
#' @param path_to_dir Path to a directory containing matrices (not the matrix folder itself). Matrix must be must in cellphylo format
#'
#' @return For each matrix in the directory, prints out stats for the number of cell type groups, the number of replicates per cell type group and identifies the cell type group with the fewest number of replicates.
#' @export
#'
#' @examples
#'
count_replicates <- function(path_to_dir){
  #setup

  dir_list <- list.dirs(path = path_to_dir, full.names=TRUE, recursive=FALSE )

  lapply(dir_list, function(x){

    #load matrix
    print("Reading in matrix...")
    mat <- Read10X(x)
    print("Matrix loaded.")
    #parse names
    species <- sapply(strsplit(colnames(mat)[1], "_"), function(v){return(v[1])})
    cluster_ids <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[3])})

    #identify unique cluster ids
    cluster_ids.unique <- unique(cluster_ids)
    n_clusters <- length(cluster_ids.unique)
    print(paste0("For ",species, ", there are ", n_clusters," unique cell type groups."))

    #count table
    print(paste0("Cell replicates tally for ",species, ":"))
    print(table(cluster_ids))

    #what is the smallest number of replicates?
    min_reps <- table(cluster_ids) %>% min()
    smallest_cluster <- names(table(cluster_ids))
    smallest_cluster <- smallest_cluster[which(table(cluster_ids) == min(table(cluster_ids)))]
    print(paste0("For ", species, ", the cell type group with the fewest replicates is ", smallest_cluster, " with ", min_reps, " cells"))

  }) #close lapply


} #close function
