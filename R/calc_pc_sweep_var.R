#' calc_pc_sweep_var
#'
#' Calculates tree length (Fig. 2A), average tip and interior edge lengths (Fig. 2B), and the star-ness score (ratio of sum of tip:interior edge lengths, Fig. 2C)
#'
#' @importFrom gtools mixedsort
#' @importFrom ape read.tree
#'
#' @param tree_dir_path Path to directory containing the trees (ie. `contml outtrees`) to calculate the PC sweep variables from.
#'
#' @return var.df A data frame containing the variables in a format suitable for plotting with cellphylo::plot_pc_sweep.
#' @export
#'
#' @examples
calc_pc_sweep_var <- function(tree_dir_path){


  tree_list <- list.files(path = tree_dir_path, full.names=TRUE, recursive=FALSE )
  tree_list <- mixedsort(tree_list)
  n_trees <- c( 1: length(tree_list))

  vars <- sapply(n_trees, function(n_trees){

    #read in tree
    tree <- read.tree(tree_list[n_trees])
    #plot(tree)

    #count number of taxa
    n.tips <- length(tree$tip.label)

    #interior edges
    #pull out interior edge lengths
    int.edges <- tree$edge.length[which(tree$edge[,1] > (n.tips) & tree$edge[,2] > (n.tips))]
    #sum of interior edge lengths
    sum.int.edges <-sum(int.edges)
    #average of interior edge lengths
    ave.int.edges <- mean(int.edges)


    #tip edges
    #pull out tip edge lengths
    tip.edges <- tree$edge.length[which(tree$edge[,1] <= (n.tips) | tree$edge[,2] <= (n.tips))]
    #sum of tip edge lengths
    sum.tip.edges <- sum(tip.edges)
    #average tip edge length
    ave.tip.edges<-mean(tip.edges)


    #tree length
    TL <- sum(tree$edge.length)


    #star measure
    #star.measure <- ave.tip.edges/ave.int.edges
    star.measure <- sum.tip.edges/sum.int.edges

    return(list(ave_int_edge = ave.int.edges, ave_tip_edge = ave.tip.edges, sum_int_edge = sum.int.edges, sum_tip_edge = sum.tip.edges, TL = TL, star_measure = star.measure))
  } #sapply function close

  ) #sapply close


  #put variables into a dataframe
  var.df <- t(vars) %>% as.data.frame()

  return(var.df)

}
