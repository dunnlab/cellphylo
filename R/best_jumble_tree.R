#' best_jumble_tree
#'
#' Identify the jumeble tree with the highest sum, mean and median of jumble scores.
#'
#' @param tree_dir_path Path to a directory containing multiple jumble trees
#' @param multiphylo A multiphylo object containing multiple jumble trees
#'
#' @return NULL Prints to consolethe index of the tree with the maximum sum, mean and median of jumble scores. To identify which tree file the index corresponds to use 'list.files(tree_dir_path)'.
#' @export
#'
#' @examples
#'
best_jumble_tree <- function(tree_dir_path, multiphylo){

  if (!missing(tree_dir_path) && missing(multiphylo)){
    #tree_files  <- list.files(tree_dir_path, full.names=TRUE)
    trees <- combine_multi_trees(tree_dir_path)
    multiphylo = combine_multi_trees(tree_dir_path)
  }

  if (missing(tree_dir_path) && !missing(multiphylo)){
    trees <- multiphylo
  }


  scored_trees <- lapply(trees, function(trees){
    focal_tree_path <- trees
    output <- score_jumble(multiphylo = multiphylo,focal_tree_obj = trees, print=FALSE )
    return(output)
  }) #function close, lapply close


  #max sum of jumble scores
  jumble_sums <- lapply(scored_trees, function(scored_trees){
    sums <- scored_trees$score_sum
  })
  index = which(jumble_sums == max(unlist(jumble_sums)))
  print(paste0("max sum of scores =  index: " , index, " score: ", max(unlist(jumble_sums))))

  #max mean of scores
  jumble_means <- lapply(scored_trees, function(scored_trees){
    mean <- scored_trees$score_mean
  })
  index = which(jumble_means == max(unlist(jumble_means)))
  print(paste0("max mean of scores =  index: " , index, " score: ", max(unlist(jumble_means))))

  #max median of scores
  jumble_medians <- lapply(scored_trees, function(scored_trees){
    median <- scored_trees$score_median
  })
  index = which(jumble_medians == max(unlist(jumble_medians)))
  print(paste0("max median of scores =  index: " , index, " score: ", max(unlist(jumble_medians))))

} #close func
