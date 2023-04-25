#' run_pca
#'
#' Centres the value of a matrix around the mean, performs PCA on this matrix and normalizes variance by the standard deviation.
#' The number of principal components to retain in the final PCA matrix may be indicated with the "n_PCs" option.
#' It is recommended that PCA be performed on the full, integrated unsubsampled matrix. The final PCA matrix may subset to a list of cells specified by the user by the `"subset" or "subset_file_path" option.
#'
#' @importFrom Seurat Read10X as.sparse
#' @importFrom stats prcomp sd
#' @importFrom utils read.table
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to matrix to perform PCA on
#' @param n_PCs Number of principal components to reduce rank to
#' @param subset A list of cell ids to subsample from the input matrix and retain in the resulting PCA matrix this function creates. If no subset list is provided all cells are kept.
#' @param subset_file_path A text file containing a single column of cell ids to subsample from the input matrix and retain in the resulting PCA matrix this function creates.
#' @param print Print out the PCA matrix
#'
#' @return mat.sparse A matrix PCA has been performed on in sparse matrix format. This matrix may have subsetted from the full PCA matrix calculated from the input matrix by number of PCs and cell subset to retain.
#' @export
#'
#' @examples
run_pca <- function(matrix_path, n_PCs=20, subset, subset_file_path, print=TRUE){

  #read in matrix
  mat <- Read10X(matrix_path)

  #orient matrix
  mat <- as.matrix(mat) %>% t()

  #pca
  pc <- prcomp(mat,center=TRUE, scale=FALSE)
  x <- pc$x
  #normalize by sd after rotation
  sds <- apply(x,2,sd)
  x_norm <- t(t(x)/sds)
  apply(x_norm, 2, sd)

  #take just the first 20 pcs
  #check n PCs < n cells and > 1 PC. Matrix was rotated so n cells = n rows
  if (n_PCs > nrow(mat)){
    stop("The number of PCs to keep must be fewer than the number of cells")
  }
  if (n_PCs < 2){
    stop("The number of PCs to keep must be greater than 1.")
  }
  x_norm <- x_norm[,1:n_PCs]
  print(paste0("Subsetted to ", ncol(x_norm), " PCs."))

  # subset x to 1 cell per cell type cluster per species
  #Was a path or object given?
  if (!missing(subset_file_path) && missing(subset)){
    selected <- read.table(subset_file_path) %>% unlist()
  }

  if (!missing(subset) && missing(subset_file_path)){
    selected <- subset %>% unlist()
  }

  if (missing(subset) && missing(subset_file_path)){
    selected <- rownames(x_norm)
  }

  x.sub <- x_norm[rownames(x_norm) %in% selected,]
  print(paste0("Subsetted ", nrow(mat), "cells to ", nrow(x.sub), " cells."))

  #write matrix
  mat.sparse <- as.sparse(x.sub)

  if (print==TRUE){
    if(!dir.exists("matrix")){
      dir.create("matrix")
    }

    #create a directory for the combined matrix
    if(!dir.exists("matrix/PCA")){
      dir.create("matrix/PCA")
    }

    #write matrix
    write10xCounts(path=paste0("contml_", nrow(x.sub), "cell_subset_", n_PCs, "PC_mtx"), x= mat.sparse, version="3")
    file.rename(paste0("contml_", nrow(x.sub), "cell_subset_", n_PCs, "PC_mtx"), paste0("matrix/PCA/contml_",nrow(x.sub), "cell_subset_", n_PCs, "PC_mtx"))

  } #close if print

  return(mat.sparse)


} #close func
