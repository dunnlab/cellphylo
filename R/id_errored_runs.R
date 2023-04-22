#' id_errored_runs
#'
#' Create a list of matrices produced by contml runs that errored out.
#' Errored contml runs lead to empty outtree files, which must be removed for downstream analyses. However, the corresponding scjackknife matrices must also be removed to preserve index number in downstream analyses.
#' This function is part of the workflow of the scjackknife analysis.
#'
#' @importFrom utils write.table

#' @param outtrees_dir Path to directly with contml outtrees.
#'
#' @return NULL prints out a list of scjackknife matrix names that produced empty tree files when run through contml.
#' @export
#'
#' @examples
id_errored_runs <- function(outtrees_dir){

  files_list <- list.files(path = outtrees_dir, full.names=TRUE, recursive=FALSE )

  empty_paths <- files_list[file.size(files_list) == 0]
  empty_files <- sapply(strsplit(empty_paths, "/"), function(v){return(v[length(v)])})

  #convert to matrices dir names
  names <- sapply(strsplit(empty_files, "_"), function(v){return(v[length(v)-1])})
  matrix_names <- paste0("scjackknife_540_final_", names,"_mtx")

  write.table(matrix_names, file = "errored_runs.txt", quote=FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)
}
