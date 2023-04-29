#' score_jumble
#'
#' Calculate the jumble score from a collection of jumble trees and plot onto a focal tree. The jumble score is the clade frequency of each split across jumble trees.
#'
#' @importFrom phangorn plotBS
#' @importFrom ape read.tree consensus
#' @importFrom stats median
#'
#' @param multiphylo A multiphylo object containing the jumble trees
#' @param focal_tree_path Path to the tree to be scored.
#' @param focal_tree_obj An object of class phylo containing the tree to be scored.
#' @param print Boolean. Print out the scored tree and 50% consensus tree
#'
#' @return output List object containing the scored tree, the 50% consensus tree, and the sum, mean and median of the scores of the focal tree
#' @export
#'
#' @examples
score_jumble <- function(multiphylo, focal_tree_path, focal_tree_obj, print){


  #multiphylo file
  jumble_trees <- multiphylo

  #tree to plot jaccknife scores onto
  #if given as a path to tree
  if (!missing(focal_tree_path) && missing(focal_tree_obj)){
    main_tree <- read.tree(focal_tree_path)
  }
  #if given as a phylo object
  if (missing(focal_tree_path) && !missing(focal_tree_obj)){
    main_tree <- focal_tree_obj
  }


  #calculate jackknife scores
  scored_tree <- plotBS(main_tree, jumble_trees, type="none", bs.col="red", method = "FBP")

  if (print==TRUE){
    #for soem reason figtree does not display tree if not midpoint-rooted. ape::write.tree still prints correctly even if unrooted. midpoint() before writing
    #print out tree
    ape::write.tree(scored_tree, file="scored_tree.tre")
  }

  #score stats
  score_sum <- sum(na.omit(scored_tree$node.label))
  score_mean <- mean(na.omit(scored_tree$node.label))
  score_median <- median(na.omit(scored_tree$node.label))

  #consensus tree
  consensus_tree <- consensus(jumble_trees,p=0.50, check.labels=TRUE)
  #plot tree fails at p=0.25
  #consensus_tree_25 <- consensus(jumble_trees,p=0.25, check.labels=TRUE)

  if (print==TRUE){
    ape::write.tree(consensus_tree, file="0.50_consensus_tree.tre")
    #write.tree(consensus_tree_25, file="jumble_92_0.25_consensus_tree.tre")
  }


  output <- list(scored = scored_tree, consensus = consensus_tree, score_sum = score_sum, score_mean = score_mean, score_median = score_median)
  return(output)
}
