#' make_pc_sweep_infiles
#'
#' @importFrom Seurat Read10X
#' @importFrom stats prcomp
#' @importFrom utils write.table
#'
#' @param matrix_path Path to matrix. For Mah and Dunn (2023) use `aqhumor_cross_species_integrated_mtx`
#'
#' @return NULL
#' @export
#'
#' @examples
make_pc_sweep_infiles <- function(matrix_path){

  #unfiltered raw data
  mat <- Read10X(matrix_path)
  mat <- as.matrix(mat) %>% t()

  #pca
  pc <- prcomp(mat, center=TRUE, scale=FALSE)
  x <- pc$x
  #normalize sd
  sds <- apply(x,2,sd)
  x.scaled <- t(t(x)/sds)
  apply(x.scaled, 2, sd)

  #c = 1 doesn't work for write_matrix_to_phylip()
  #col_range <- 2:ncol(x)
  #make only 100 trees
  col_range <- 2:100

  lapply(col_range, function(c){
    x.sweep.matrix <- x.scaled[,1:c] %>% as.matrix()
    #format to phylip
    phylip <- make_infile(matrix = x.sweep.matrix, print=FALSE)
    #write table
    write.table(phylip, file=paste0("infile_pca_var_norm_pc_sweep_",c), quote=FALSE, sep=" ",   row.names=FALSE, col.names=FALSE, na= "")
  } #close function
  ) #close lapply

}
