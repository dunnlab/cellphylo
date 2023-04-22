#' color_tree
#' 
#' @import ggtree 
#' @importFrom phangorn midpoint
#' @importFrom ape read.tree
#' @importFrom tibble column_to_rownames
#' @importFrom grDevices pdf dev.off
#' 
#' @param tree_path Path to tree
#' @param print Print out the colored trees as pdf
#'
#' @return A list containing two trees: 'tip_tree' = the tree with tip labels colored by cell type and 'bar_tree': the tree with cell types indicated with a colored bar
#' @export
#'
#' @examples
color_tree = function(tree_path, print){

  tree <- read.tree(tree_path)
  #midpoint root tree
  tree <- midpoint(tree)
  #parse tips to make annotation df
  #designed for scjack_outtrees_labels
  celltype_ann <- sapply(strsplit(tree$tip.label, "_"), function(v){return(v[4])})
  #annotation df to attach to ggtree object = tips = column 1
  annotation.df.1 <- data.frame("tip" = tree$tip.label, "annotation"= celltype_ann)
  #annotation df to use in gheatmap - tips = rownames
  annotation.df.2 <- column_to_rownames(annotation.df.1, var="tip")
  
  
  #colour tips by annotation
  tree.tip_cols <- ggtree(tree, aes(color=annotation)) %<+% annotation.df.1 + geom_tiplab() 
  
  #label tip by colour bar
  tree.gg <- ggtree(tree)  %<+% annotation.df.2
  tree.heatmap <- gheatmap(tree.gg, annotation.df.2,  width=0.01, colnames=FALSE, legend="cell type", offset=1.5) + geom_tiplab()
  
  if (print==TRUE){
    pdf(file="color_tip.pdf", width=20, height=20)
    print(tree.tip_cols)
    dev.off()
    
    pdf(file="color_bar.pdf", width=20, height=20)
    print(tree.heatmap)
    dev.off()
    
  }#close print if
  
  return(list(tip_tree = tree.tip_cols, bar_tree = tree.heatmap))
} #close function
