#' label_scjack_trees_54
#' 
#' Comparison and calculation of clade frequencies across scjackknife trees require that the tips share identical names
#' This script labels scjackknife outtrees by 1) their full tip names and 2) a common label - a tip name that contains only species and cell type cluster group id that can be matched acros scjackknife trees
#' 
#' @importFrom ape write.tree
#' 
#' @param n_trees The number of trees to be relabelled. 
#' @param matrix_dir_path The path to the scjackknife matrices that were used to create the scjackknife infiles. This directory must have matrices that produced empty outtree files (ie contml run errored) removed using cellphylo::id_errored_runs.
#' @param tree_dir_path The path to the scjakknife outtrees produced by contml. This must have empty outtree files removed (see scjackknife tutorial).
#'
#' @return output A list containing scjackknife trees with full tip name labels ("labelled_tree") and scjackknife trees with comomn labels ("common_labels_tree)
#' @export
#'
#' @examples

label_scjack_trees_54<- function(n_trees, matrix_dir_path, tree_dir_path){

  tree_list <- list.files(path = tree_dir_path, full.names=TRUE, recursive=FALSE )  
  matrix_list <- list.dirs(path = matrix_dir_path, full.names=TRUE, recursive=FALSE)
  
  sapply(n_trees, function(tree_index){
    labelled_tree <- label_tree(matrix_path = matrix_list[tree_index], tree_path = tree_list[tree_index], print=FALSE)
    
    #parse out tree name
    file_name <- sapply(strsplit(tree_list[tree_index], "/"), function(v){return(v[length(v)])})
    index_num <- sapply(strsplit(file_name, "_"), function(v){return(v[length(v)-1])})
    suffix <- sapply(strsplit(file_name, "_"), function(v){return(v[length(v)])})
    
    #print out trees with human readable tips
    write.tree(labelled_tree, file=paste0("scjackknife_540_outtrees_labels_", index_num, "_", suffix))
    
    #bootstrapping requires common label across jackknife trees so that tips correspond
    species <- sapply(strsplit(labelled_tree$tip.label, "_"), function(v){return(v[1])})
    cluster_id <- sapply(strsplit(labelled_tree$tip.label, "_"), function(v){return(v[3])})
    new_label <- paste0(species, "_", cluster_id)
    common_labels_tree <- labelled_tree
    common_labels_tree$tip.label <- new_label
    #print out common labels tree
    write.tree(common_labels_tree, file=paste0("scjackknife_540_outtrees_common_labels_", index_num, "_", suffix))
    
    output <- list("labelled_tree" = labelled_tree, "common_labels_tree" = common_labels_tree)
    return(output)
    
  } #close sapply function
  
  
  )#close sapply
  
} 