#' label_scjack_trees
#'
#' Comparison and calculation of TBE across scjackknife trees require that the tips share identical names
#' This script labels scjackknife outtrees by 1) their full tip names and 2) a common label - a tip name that contains only species and cell type cluster group id that can be matched acros scjackknife trees
#' Previously named label_scjack_trees_54.
#'
#' @importFrom ape write.tree
#'
#' @param n_trees The number of trees to be relabelled.
#' @param matrix_dir_path The path to the scjackknife matrices that were used to create the scjackknife infiles. This directory must have matrices that produced empty outtree files removed using cellphylo::id_errored_runs.
#' @param tree_dir_path The path to the scjakknife outtrees produced by contml. This must have empty outtree files removed.
#' @param print Boolean. Print out trees with 1) full labels and 2) common labels
#' @param file_name_prefix The prefix for the outtree files to produce.
#' @return output A list containing scjackknife trees with full tip name labels ("labelled_tree") and scjackknife trees with comomn labels ("common_labels_tree)
#' @export
#'
#' @examples

label_scjack_trees<- function(n_trees, matrix_dir_path, tree_dir_path, print=TRUE, file_name_prefix){
   if (print==TRUE){
     if(!dir.exists(paste0(file_name_prefix, "_label"))){
       dir.create(paste0(file_name_prefix, "_label"))
     }
     if(!dir.exists(paste0(file_name_prefix, "_common_label"))){
       dir.create(paste0(file_name_prefix, "_common_label"))
     }
   }

  tree_list <- list.files(path = tree_dir_path, full.names=TRUE, recursive=FALSE )
  matrix_list <- list.dirs(path = matrix_dir_path, full.names=TRUE, recursive=FALSE)



  lapply(c(1:n_trees), function(tree_index){
    #label the tree with the full cell ids
    labelled_tree <- label_tree(matrix_path = matrix_list[tree_index], tree_path = tree_list[tree_index], print=FALSE)

    #parse out cell id
    file_name <- sapply(strsplit(tree_list[tree_index], "/"), function(v){return(v[length(v)])})
    index_num <- sapply(strsplit(file_name, "_"), function(v){return(v[length(v)-1])})
    suffix <- sapply(strsplit(file_name, "_"), function(v){return(v[length(v)])})

    if (print==TRUE){
    #print out trees with human readable tips. Full cell ids.
    write.tree(labelled_tree, file=paste0(file_name_prefix,"_labels_", index_num, "_", suffix))
    file.rename(paste0(file_name_prefix,"_labels_", index_num, "_", suffix), paste0(file_name_prefix, "_label/",file_name_prefix,"_labels_", index_num, "_", suffix))
    }

    #bootstrapping requires common label across jackknife trees so that tips correspond
    #make common labels
    species <- sapply(strsplit(labelled_tree$tip.label, "_"), function(v){return(v[1])})
    cluster_id <- sapply(strsplit(labelled_tree$tip.label, "_"), function(v){return(v[3])})
    new_label <- paste0(species, "_", cluster_id)
    common_labels_tree <- labelled_tree
    common_labels_tree$tip.label <- new_label

    if (print==TRUE){
    #print out common labels tree
    write.tree(common_labels_tree, file=paste0(file_name_prefix, "_common_labels_", index_num, "_", suffix))
    file.rename(paste0(file_name_prefix, "_common_labels_", index_num, "_", suffix), paste0(file_name_prefix, "_common_label/",file_name_prefix, "_common_labels_", index_num, "_", suffix))
    }

    output <- list("labelled_tree" = labelled_tree, "common_labels_tree" = common_labels_tree)

    return(output)

  } #close apply function

  )#close apply

}
