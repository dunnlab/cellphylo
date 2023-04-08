#' label_tree
#'
#' Replace the taxon number labels in the contml outtree with the cell ids in the input matrix used to make the contml infile.
#'
#' @importFrom Seurat Read10X
#' @importFrom ape read.tree write.tree
#' @importFrom dplyr left_join
#'
#' @param matrix_path Path to matrix from which the contml infile was created. The taxon numbers of the infile correspond to the matrix cell ids in sequence.
#' @param tree_path Path to the raw phylogeny produced by contml ("outtree" file)
#' @param file_name Name of tree file of labelled tree
#' @param print Print out the tree file of the labelled tree
#'
#' @return A phylo object of the labelled tree
#' @export
#'
#' @examples
label_tree = function(matrix_path,tree_path, file_name, print){

  mat <- Read10X(matrix_path)

  #connect cell ids to index
  key.df <- data.frame(cell_id = rownames(mat), index = as.character(c(1:nrow(mat))))

  #load tree
  tree <- read.tree(tree_path)
  plot(tree)
  #find tip label
  tips.df <- data.frame(index = tree$tip.label)

  #connect tip id to cell id
  new_tip.df <- left_join(tips.df, key.df, by="index")

  #replace old tip indices with new tip cell id labels
  tree$tip.label <- new_tip.df$cell_id
  plot(tree)

  if (print==TRUE){
    write.tree(tree, file=file_name)
  }

  return(tree)
}
