#' make_scjack_sample
#'
#' Create resampled matrices that contain 1 replicate cell per cell type group from a larger matrix that contains multiple replicates. This is used for running the scjackknife analysis.
#' Previously named make_scjack_sample54
#'
#' @importFrom utils write.table
#' @importFrom DropletUtils write10xCounts
#' @importFrom Seurat Read10X
#'
#' @param n_samples Number of scjackknife samples to create.
#' @param subsample_matrix_path Path to the matrix to sample from.
#' @param full_matrix_path Path to full unsubsetted matrix. PCA should already have been performed and PCs subsetted.
#' @return pca_matrix List of resampled single jackknife matrices.
#' @export
#'
#' @examples
make_scjack_sample <- function(n_samples, subsample_matrix_path, full_matrix_path, file_name){
   #load the full 919 cell 20 PC matrix for subsetting down below
   full_pca_mat <- Read10X(full_matrix_path)

   #make directories to organize
   if(!dir.exists("scjackknife")){
     dir.create("scjackknife")
   }

   if(!dir.exists(paste0("scjackknife/input_", file_name))){
     dir.create(paste0("scjackknife/input_", file_name))
   }

   if(!dir.exists(paste0("scjackknife/input_", file_name,"/scjack_selected_cells"))){
     dir.create(paste0("scjackknife/input_", file_name,"/scjack_selected_cells"))
   }

   if(!dir.exists(paste0("scjackknife/input_", file_name,"/scjack_matrices"))){
     dir.create(paste0("scjackknife/input_", file_name,"/scjack_matrices"))
   }

   if(!dir.exists(paste0("scjackknife/input_", file_name, "/scjack_infiles"))){
     dir.create(paste0("scjackknife/input_", file_name, "/scjack_infiles"))
   }

  sapply(n_samples, function(n_samples){

    #draw random subset of cells. 1 cell per cell type cluster per speciesu
    sampled_cells <- subset_combined_matrix(matrix_path = subsample_matrix_path, sample_size=1, print=FALSE)
    #record list of selected cells
    subset_list <- colnames(sampled_cells)

    #rank reduce matrix
    #pca_matrix <- run_pca(matrix_path = matrix_path, n_PCs = 20, subset = subset_list, print=FALSE)
    pca_matrix <- full_pca_mat[rownames(full_pca_mat) %in% subset_list,]

    #make infile
    phylip <- make_infile(matrix = pca_matrix, print = FALSE)

    ## print out scjackknife input files
    ## print outside of the two above functions so can have sequential file names

    ## record the subsampled cell ids
    write.table(subset_list, paste0("scjackknife_", file_name,"_selected_cells_",n_samples, ".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    file.rename(paste0("scjackknife_", file_name,"_selected_cells_",n_samples, ".txt"), paste0("scjackknife/input_", file_name,"/scjack_selected_cells/scjackknife_", file_name,"_selected_cells_",n_samples, ".txt"))
    ## write out subasmpled matrix
    write10xCounts(paste0("scjackknife_", file_name,"_",n_samples, "_mtx"),  pca_matrix, version="3")
    file.rename(paste0("scjackknife_", file_name,"_",n_samples, "_mtx"), paste0("scjackknife/input_", file_name,"/scjack_matrices/scjackknife_", file_name,"_",n_samples, "_mtx"))
    ## phylip infile
    write.table(phylip, file=paste0("infile_scjackknife_", file_name,"_", n_samples), quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE, na= "")
    file.rename(paste0("infile_scjackknife_", file_name,"_", n_samples), paste0("scjackknife/input_", file_name,"/scjack_infiles/infile_scjackknife_", file_name,"_", n_samples))

    return(pca_matrix)

  } #close function

  )#close sapply

} #close scjackknife_make_input()
