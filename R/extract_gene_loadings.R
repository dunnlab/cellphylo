#' extract_gene_loadings
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr left_join select 
#' @importFrom Seurat Read10X
#' @importFrom stats prcomp 
#' @importFrom tibble rownames_to_column 
#' 
#' @param ann_table Annotation table connecting gene ids to human gene symbols.
#' @param matrix_path Path to matrix in 10x format to do PCA on and extract gene loadings. 
#'
#' @return gene_loadings_output A list containing 1) 'by_most_Positive': a data frame with gene symbols have been ordered by most positive loading for each PC 2) "by_most_Negative": A data frame with gene symbol have been ordered by mmost negative loading for each PC, "PCA_output": the object produced by performing PCA on the input matrix with prcomp.
#' @export
#'
#' @examples
extract_gene_loadings <- function(ann_table, matrix_path){
  
#load gene id to gene symbol ann file
ann_table <- read_tsv(ann_table)
#change underscores to dashes to match matrix
ann_table$gene_identifier <- gsub("_", "-", ann_table$gene_identifier)
#just use the human gene symbols
ann_table <- select(ann_table, c("gene_identifier", "human"))
colnames(ann_table) <- c("gene_identifier", "gene_symbol")

#load matrix
mat <- Read10X(matrix_path)

#replace gene ids with gene symbols
mat_names <- rownames(mat) %>% as.data.frame()
colnames(mat_names) <- "gene_identifier"
new_mat_names <- left_join(mat_names, ann_table, by = "gene_identifier")
new_mat_names <- new_mat_names %>% select("gene_symbol")
rownames(mat) <- new_mat_names[["gene_symbol"]]

#PCA
#must transpose for prcomp
mat <- as.matrix(mat) %>% t()
#prcomp
pca <- prcomp(mat,center=TRUE, scale=FALSE)
#pca$rotation are the gene loadings for each pc
loadings <- pca$rotation %>% as.data.frame()

#prep loadings df for lapply fun
loadings <- rownames_to_column(loadings, var="gene_symbol")
#index for lapply fun. First col of loadings is not a PC
n_pc <- 1:(ncol(loadings)-1)


#order each PC by most positive and most negative gene loadings
order  <- lapply(n_pc, function(pc_index){
  #name of PC column are we working iwth
  PCn <- paste0("PC", pc_index)
  #isolate gene symbol col and PC col
  pc_column <- loadings[,c("gene_symbol", PCn)]
  # "-" sorts by descending (highest value first)
  by_most_pos <- pc_column[order(-pc_column[,2]),] %>% select(gene_symbol)
  colnames(by_most_pos) <- PCn
  by_most_neg <-  pc_column[order(pc_column[,2]),] %>% select(gene_symbol)
  colnames(by_most_neg) <- PCn
  
  output <- list("positive" = by_most_pos, "negative" = by_most_neg)
  
  return(output)
})
positive <- sapply(order, "[[", "positive")
positive_df <- do.call(cbind, positive)

negative <- sapply(order, "[[", "negative")
negative_df <- do.call(cbind, negative)

gene_loadings_output <- list("by_most_Positive" = positive_df, "by_most_Negative" = negative_df, "PCA_output" = pca)

return(gene_loadings_output)

}