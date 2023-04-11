#' combine_multi_trees
#'
#' @importFrom TreeTools as.multiPhylo
#' @importFrom ape read.tree
#' @param tree_dir_path Path to directory containing multiple trees to combine
#'
#' @return multi_tree A multiphylo object of all trees in the tree_dir_path directory
#' @export
#'
#' @examples
#'
combine_multi_trees <- function(tree_dir_path){

  tree_paths <- list.files(tree_dir_path, full.names=TRUE)

  multi_tree <- lapply(tree_paths, function(tree_paths){
    tree <- read.tree(tree_paths)
    return(tree)
  })

  #convert multi_tree list to multiphylo object
  #as.multiPhylo not to be mixed up with phytools::as.multiPhylo
  multi_tree <- TreeTools::as.multiPhylo(multi_tree)

  return(multi_tree)
}
