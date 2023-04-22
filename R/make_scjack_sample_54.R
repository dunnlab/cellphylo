#' make_scjack_sample_54
#' Create resampled matrices that contain 1 replicate cell per cell type group from a larger matrix that contains multiple replicates. This is used for running the scjackknife analysis.
#'
#' @importFrom utils write.table
#' @importFrom DropletUtils write10xCounts
#' @param n_samples Number of scjackknife samples to create.
#' @param matrix_path Path to the matrix to sample from.
#'
#' @return pca_matrix List of resampled single jackknife matrices.
#' @export
#'
#' @examples
make_scjack_sample_54 <- function(n_samples, matrix_path){
  sapply(n_samples, function(n_samples){

    #draw random subset of cells. 1 cell per cell type cluster per speciesu
    sampled_cells <- subset_combined_matrix(matrix_path = matrix_path, sample_size=1, print=FALSE)
    #record list of selected cells
    subset_list <- colnames(sampled_cells)

    #rank reduce matrix
    pca_matrix <- run_pca(matrix_path = matrix_path, n_PCs = 20, subset = subset_list, print=FALSE)

    ## print out scjackknife input files
    ## print outside of the two above functions so can have sequential file names

    ## record the subsampled cell ids
    write.table(subset_list, paste0("scjackknife_540_final_selected_cells_",n_samples, ".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    ## write out subasmpled matrix
    write10xCounts(paste0("scjackknife_540_final_",n_samples, "_mtx"),  pca_matrix, version="3")
    ## phylip infile
    phylip <- make_infile(matrix = pca_matrix, print = FALSE)
    write.table(phylip, file=paste0("infile_scjackknife_540_final_", n_samples), quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE, na= "")


    return(pca_matrix)

  } #close function

  )#close sapply

} #close scjackknife_make_input()
