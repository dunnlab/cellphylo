#' best_scjack_tree
#'
#' @param tree_dir_path Path to a directory containing multiple jumble trees
#' @param multiphylo A multiphylo object containing multiple jumble trees
#'
#' @return NULL Prints out the index of the tree with the maximum sum, mean and median of jumble scores. To identify which tree file the index corresponds to use 'list.files(tree_dir_path)'.
#' @export
#'
#' @examples
best_scjack_tree <- function(tree_dir_path, multiphylo){

  if (!missing(tree_dir_path) && missing(multiphylo)){
    #tree_files  <- list.files(tree_dir_path, full.names=TRUE)
    trees <- combine_multi_trees(tree_dir_path)
    scjack_multiphylo = combine_multi_trees(tree_dir_path)
  }

  if (missing(tree_dir_path) && !missing(multiphylo)){
    trees <- multiphylo
    scjack_multiphylo <- multiphylo
  }


  scored_trees <- lapply(trees, function(trees){
    focal_tree_path <- trees
    output <- score_scjack(scjack_multiphylo = scjack_multiphylo,focal_tree_obj = trees, print=FALSE )
    return(output)
  }) #function close, lapply close


  #max sum of scjack scores
  scjack_sums <- lapply(scored_trees, function(scored_trees){
    sums <- scored_trees$score_sum
  })
  index = which(scjack_sums == max(unlist(scjack_sums)))
  print(paste0("max sum of scores =  index: " , index, " score: ", max(unlist(scjack_sums))))

  #max mean of scores
  scjack_means <- lapply(scored_trees, function(scored_trees){
    mean <- scored_trees$score_mean
  })
  index = which(scjack_means == max(unlist(scjack_means)))
  print(paste0("max mean of scores =  index: " , index, " score: ", max(unlist(scjack_means))))

  #max median of scores
  scjack_medians <- lapply(scored_trees, function(scored_trees){
    median <- scored_trees$score_median
  })
  index = which(scjack_medians == max(unlist(scjack_medians)))
  print(paste0("max median of scores =  index: " , index, " score: ", max(unlist(scjack_medians))))

} #close func
