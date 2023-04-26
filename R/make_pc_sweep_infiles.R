#' make_pc_sweep_infiles
#'
#' This function produces the input files for the PC sweep experiment in Mah & Dunn (2023) (Fig.2).
#' It takes in the 919 cell 919 PC matrix (Fig. S1, step 7) and subsets the matrix to an increasing number of PCs. These matrix subsets are printed out as `contml` infiles.
#' If `calculate_PCA = TRUE` the function will first perform PCA on the matrix before subsetting and creating the infiles.
#'
#' @param matrix_path Path to a matrix that has had PCA performed on it.
#' @param run_PCA Boolean. If TRUE, calculate PCA on the input matrix before subsetting and producing infiles (i.e. the input matrix has not had PCA performed on it yet)
#' @param PC_range A range indicating the number of PCs to subset from the input matrix.
#'
#' @return This function prints out the `contml` infiles for the PC sweep experiment
#' @export
#'
#' @examples
#'
make_pc_sweep_infiles <- function(matrix_path, run_PCA = FALSE, PC_range=c(3:100)){


  mat <- Read10X(matrix_path) %>% as.matrix()

  if(run_PCA == TRUE){
    mat <- t(mat)

    #pca
    pc <- prcomp(mat, center=TRUE, scale=FALSE)
    x <- pc$x
    #normalize sd
    sds <- apply(x,2,sd)
    x.scaled <- t(t(x)/sds)
    apply(x.scaled, 2, sd)
    mat <- x.scaled
  }
  #set up directories to deposit infiles into
  if(!dir.exists("PC_sweep")){
    dir.create("PC_sweep")
  }

  #create a directory for the combined matrix
  if(!dir.exists("PC_sweep/infiles")){
    dir.create("PC_sweep/infiles")
  }


  lapply(PC_range, function(c){
    x.sweep.matrix <- mat[,1:c]
    #format to phylip
    phylip <- make_infile(matrix = x.sweep.matrix, print=FALSE)
    #write table
    write.table(phylip, file=paste0("PC_sweep/infiles/infile_pc_sweep_",c), quote=FALSE, sep=" ",   row.names=FALSE, col.names=FALSE, na= "")
    print(paste0("infile ", c))
  } #close function
  ) #close lapply

}
