#' combine_matrices
#'
#' Takes in multiple species matrices and combines them into a single multi-species matrix combined by gene id. This combined matrix is NOT cross-species normalized yet.
#'
#' @import dplyr
#' @importFrom Seurat Read10X as.sparse
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom DropletUtils write10xCounts
#' @importFrom stats na.omit
#'
#' @param matrix_list  A list of the full paths to the matrices you want to combine. Cell ids must be in cellphylo format.
#' @param ortholog_map_path Path to wrangled Ensembl orthologs file. The custom script for wrangling the raw Ensembl orthologs file can be found in wrangle_Ensembl_orthologs.Rmd.
#' @param combined_mat_name Name for the combined matrix this function produces. Cell ids are in cellphylo format.
#'
#' @return full_join_sparse The cross-species combined matrix in sparse matrix format
#' @export
#'
#' @examples
combine_matrices <- function(matrix_list, ortholog_map_path, combined_mat_name){

  ensembl_map_final <- read_tsv(ortholog_map_path)
  species_cols <- ensembl_map_final[, !names(ensembl_map_final) %in% "gene_identifier"] %>% colnames()
  gene_ont <- pivot_longer(ensembl_map_final, cols =  all_of(species_cols))
  colnames(gene_ont) <- c("gene_identifier", "species", "gene")

  #initialize an empty list
  list_matrices <- list()

  for (i in 1:length(matrix_list)){
    #load matrix
    mat <- Read10X(matrix_list[i])
    mat <- as.data.frame(mat)
    mat <- rownames_to_column(mat, var = "gene")

    #get species name
    current_species <- sapply(strsplit(colnames(mat)[2], "_"), function(v){return(v[1])})


    #filter gene_ont for the current species matrix
    gene_ont_sub <- filter(gene_ont, species == current_species)
    gene_ont_sub <- gene_ont_sub %>% select(gene_identifier, gene)

    #inner join - add gene_identifiers to mat
    joined.matrix <- inner_join(mat, gene_ont_sub, by="gene")
    #remove column with gene symbols - will use gene_identifiers instead
    joined.matrix <- subset(joined.matrix, select = -gene)

    if (i == 1) { full_join <- joined.matrix }
    if (i > 1) { full_join <- full_join(full_join, joined.matrix) }

    list_matrices[[i]]  <- joined.matrix


  } #end forloop


  #format matrix
  full_join <- column_to_rownames(full_join, var = "gene_identifier")
  #remove NAs
  full_join_noNA <- na.omit(full_join)

  full_join_sparse <- as.sparse(full_join_noNA)
  write10xCounts(combined_mat_name, full_join_sparse, version="3")

  return(full_join_sparse)

} #end function
