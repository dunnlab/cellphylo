
#' make_infile
#'
#' Converts a matrix to contml (PHYLIP) infile format.
#'
#' @importFrom Seurat Read10X
#' @importFrom utils write.table

#' @param matrix Matrix to create contml infile from. Use this if you have the matrix as an R object.
#' @param matrix_path Path to matrix to create contml infile from. Use this if you have a printed matrix in 10x feature-barcode format.
#' @param print Print the infile
#'
#' @return An object representing a matrix formatted in contml infile format. Contml requires a printed infile so print=TRUE is recommended.
#' @export
#'
#' @examples
make_infile = function(matrix, matrix_path, print){

  if (!missing(matrix_path) && missing(matrix)){
    mat <- Read10X(matrix_path)
  }

  if (!missing(matrix) && missing(matrix_path)){
    mat <- matrix
  }

  ntaxa <- nrow(mat)
  nchar <- ncol(mat)
  NA_col <- rep(NA, times=nrow(mat))
  phylip=""
  phylip <- cbind(1:nrow(mat), NA_col, NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col, mat)
  phylip_final <- rbind(c(rep(NA, times=4), ntaxa, rep(NA, times=3), nchar, rep(NA, times=ncol(phylip)-9)), phylip)

  if (print == TRUE){
    print(write.table(phylip_final, file="infile", quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE, na= ""))}

  return(phylip_final)

}
