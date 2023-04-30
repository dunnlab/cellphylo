#' id_errored_runs
#'
#' When `contml` runs error, it produces an empty `outtree` and `outfile`. We remove these, but now the index between the `contml` output files and scjack matrices do not match, leading to downstream analysis issues.
#' This function makes a list of matrices that produced an error when `contml` was run. We will then use this list to remove these matrices and proceed with downstream analyses.
#' This function is part of the workflow of the scjackknife analysis.
#'
#'
#' @importFrom utils write.table

#' @param outtrees_dir Path to directory with contml `outtrees`. This directory contains all `outtree` files created by contml, including for runs that errored.
#' @param matrix_name_pattern Prefix for matrix name. This function assumes that your matrix names end in `index_num_mtx` eg. for `scjackknife_540_TBE_1_mtx`, the matrix_name_pattern is `scjackknife_540_TBE`
#'
#' @return prints out a list of scjackknife matrix names that produced empty tree files when run through contml.
#' @export
#'
#' @examples
#'
id_errored_runs <- function(outtrees_dir, matrix_name_pattern){

  files_list <- list.files(path = outtrees_dir, full.names=TRUE, recursive=FALSE )

  empty_paths <- files_list[file.size(files_list) == 0]
  empty_files <- sapply(strsplit(empty_paths, "/"), function(v){return(v[length(v)])})

  #convert to matrices dir names
  names <- sapply(strsplit(empty_files, "_"), function(v){return(v[length(v)-1])})
  matrix_names <- paste0(matrix_name_pattern,"_", names,"_mtx")

  write.table(matrix_names, file = "errored_runs.txt", quote=FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)
}
