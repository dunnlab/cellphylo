#' combine_matrices
#'
#' Takes in multiple species' matrices and combines them into a single multi-species matrix combined by gene orthology. This combined matrix is NOT cross-species normalized yet.
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
#' @param ortholog_map_path Path to wrangled Ensembl orthologs file.
#' @param print Boolean. Print out the combined matrix.
#'
#' @return full_join_sparse The combined matrix in sparse matrix format. If print=TRUE, it will print out the combined matrix and deposit it in `matrix/cross-species_analysis/combined_by_orthology`
#' @export
#'
#' @examples
#'
combine_matrices <- function(matrix_list, ortholog_map_path, print=TRUE){

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

  if(print==TRUE){

    #create a matrix directory if doesn't already exist
    if(!dir.exists("matrix")){
      dir.create("matrix")
    }

    #create a directory for the combined matrix
    if(!dir.exists("matrix/combined_by_orthology")){
      dir.create("matrix/combined_by_orthology")
    }

  write10xCounts("aqhumor_combined_mtx", full_join_sparse, version="3")
  #move matrix to its own folder
  file.rename("aqhumor_combined_mtx", "matrix/combined_by_orthology/aqhumor_combined_mtx")

  }
  return(full_join_sparse)

} #end function
