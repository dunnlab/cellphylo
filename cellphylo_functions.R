library(magrittr)
library(dplyr)
library(utils)
library(DropletUtils)
library(R.utils)
library(Seurat)
library(tibble)
library(readr)
library(tidyr)
library(stats)
library(factoextra)
library(ggplot2)
library(ape)
library(TreeTools)
library(gtools)
library(phangorn)
library(ggtree)


#' best_jumble_tree
#'
#' Identify the jumeble tree with the highest sum, mean and median of jumble scores.
#'
#' @param tree_dir_path Path to a directory containing multiple jumble trees
#' @param multiphylo A multiphylo object containing multiple jumble trees
#'
#' @return NULL Prints to consolethe index of the tree with the maximum sum, mean and median of jumble scores. To identify which tree file the index corresponds to use 'list.files(tree_dir_path)'.
#' @export
#'
#' @examples
#'
best_jumble_tree <- function(tree_dir_path, multiphylo){
  
  if (!missing(tree_dir_path) && missing(multiphylo)){
    #tree_files  <- list.files(tree_dir_path, full.names=TRUE)
    trees <- combine_multi_trees(tree_dir_path)
    multiphylo = combine_multi_trees(tree_dir_path)
  }
  
  if (missing(tree_dir_path) && !missing(multiphylo)){
    trees <- multiphylo
  }
  
  
  scored_trees <- lapply(trees, function(trees){
    focal_tree_path <- trees
    output <- score_jumble(multiphylo = multiphylo,focal_tree_obj = trees, print=FALSE )
    return(output)
  }) #function close, lapply close
  
  
  #max sum of jumble scores
  jumble_sums <- lapply(scored_trees, function(scored_trees){
    sums <- scored_trees$score_sum
  })
  index = which(jumble_sums == max(unlist(jumble_sums)))
  print(paste0("max sum of scores =  index: " , index, " score: ", max(unlist(jumble_sums))))
  
  #max mean of scores
  #  jumble_means <- lapply(scored_trees, function(scored_trees){
  #    mean <- scored_trees$score_mean
  #  })
  #  index = which(jumble_means == max(unlist(jumble_means)))
  #  print(paste0("max mean of scores =  index: " , index, " score: ", max(unlist(jumble_means))))
  
  #max median of scores
  # jumble_medians <- lapply(scored_trees, function(scored_trees){
  #    median <- scored_trees$score_median
  #  })
  #  index = which(jumble_medians == max(unlist(jumble_medians)))
  #  print(paste0("max median of scores =  index: " , index, " score: ", max(unlist(jumble_medians))))
  
} #close func

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

#' combine_multi_trees
#'
#' @importFrom TreeTools as.multiPhylo
#' @importFrom ape read.tree
#' @param tree_dir_path Path to directory containing the trees to combine into a multiphylo object. These tree files should be in newick format.
#'
#' @return multi_tree A multiphylo object of all trees in the tree_dir_path directory
#' @export
#'
#' @examples
#'
combine_multi_trees <- function(tree_dir_path){
  
  tree_paths <- list.files(tree_dir_path, full.names=TRUE)
  
  multi_tree <- lapply(tree_paths, function(tree_paths){
    tree <- read.tree(tree_paths)
    return(tree)
  })
  
  #convert multi_tree list to multiphylo object
  #as.multiPhylo not to be mixed up with phytools::as.multiPhylo
  multi_tree <- TreeTools::as.multiPhylo(multi_tree)
  
  return(multi_tree)
}

#' convert_to_10x
#'
#' Converts a UMI count matrix in CSV format to 10x-style feature-barcode matrices.
#'
#' @importFrom utils read.csv
#' @importFrom R.utils gunzip
#' @importFrom tibble column_to_rownames
#' @importFrom DropletUtils write10xCounts
#' @importFrom Seurat as.sparse
#'
#' @param path_to_csv Path to csv file
#' @param species_name Name of species matrix belongs to.
#'
#' @return sparse The matrix in sparse matrix format and a directory with a matrix in 10x feature-barcode format (version 3)
#' @export
#'
#' @examples
#'
convert_to_10x <- function(path_to_csv, species_name){
  
  #gunzip
  unzipped <- gunzip(path_to_csv, remove=FALSE)
  #read csv
  csv <- read.csv(unzipped)
  #make rownames gene ids
  csv <- column_to_rownames(csv, var = "X")
  #convert to sparse
  sparse <- as.sparse(csv)
  #write new 10x matrix
  write10xCounts(paste0("matrix_10x_", species_name), sparse, version="3")
  #create a folder for matrices, if not already present
  if(!dir.exists("matrix")){
    dir.create("matrix")
  }
  #create the 10x folder if not already present
  if(!dir.exists("matrix/10x")){
    dir.create("matrix/10x")
  }
  #move this file into the 10x folder
  file.rename(paste0("matrix_10x_", species_name), paste0("matrix/10x/matrix_10x_", species_name))
  #remove expanded file
  remove_file <- gsub(".gz", "", path_to_csv)
  file.remove(remove_file)
  
  return(sparse)
} #close fun


#' count_replicates
#'
#' This function counts the number of replicates per cell type group.
#'
#' @importFrom Seurat Read10X
#'
#'
#' @param path_to_dir Path to a directory containing matrices (not the matrix folder itself). Matrix must be must in cellphylo format
#'
#' @return For each matrix in the directory, prints out stats for the number of cell type groups, the number of replicates per cell type group and identifies the cell type group with the fewest number of replicates.
#' @export
#'
#' @examples
#'
count_replicates <- function(path_to_dir){
  #setup
  
  dir_list <- list.dirs(path = path_to_dir, full.names=TRUE, recursive=FALSE )
  
  lapply(dir_list, function(x){
    
    #load matrix
    print("Reading in matrix...")
    mat <- Read10X(x)
    print("Matrix loaded.")
    #parse names
    species <- sapply(strsplit(colnames(mat)[1], "_"), function(v){return(v[1])})
    cluster_ids <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[3])})
    
    #identify unique cluster ids
    cluster_ids.unique <- unique(cluster_ids)
    n_clusters <- length(cluster_ids.unique)
    print(paste0("For ",species, ", there are ", n_clusters," unique cell type groups."))
    
    #count table
    print(paste0("Cell replicates tally for ",species, ":"))
    print(table(cluster_ids))
    
    #what is the smallest number of replicates?
    min_reps <- table(cluster_ids) %>% min()
    smallest_cluster <- names(table(cluster_ids))
    smallest_cluster <- smallest_cluster[which(table(cluster_ids) == min(table(cluster_ids)))]
    print(paste0("For ", species, ", the cell type group with the fewest replicates is ", smallest_cluster, " with ", min_reps, " cells"))
    
  }) #close lapply
  
  
} #close function

#' cross_species_integration
#'
#' Integrate counts across species using the Seurat workflow. The default values for this function are the parameters for the Seurat workflow performed by Mah & Dunn (2023).
#'
#' @import Seurat
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to combined matrix. This combined matrix features cells from multiple species but is not yet cross-species integrated. This matrix was created using `cellphylo::combine_matrices`. Cell ids must be in cellphylo format.
#' @param min_cells input for the min.cells option of Seurat::CreateSeuratObject. Default is 3.
#' @param min_features input for the min.features option of Seurat::CreateSeuratObject. Default is 200.
#' @param k_weight The k.weight to use with Seurat::IntegrateData. This must be less than or equal to the number of cells in the smallest batch. Default is 100.
#' @param n_pcs input for the npcs option of Seurat::RunPCA. Default is 30.
#' @param UMAP_dim input for the dim option of Seurat:: RunUMAP. Default is 1:30.
#' @param FindNeighbors_dim input for the dim option of Seurat::FindNeighbors. Default is 1:30.
#' @param FindClusters_res input for the res option of Seurat::FindClusters. Default is 0.5.
#' @param print Boolean. Print out the matrix and Seurat object of cross-species integrated data.
#'
#' @return combined The Seurat object containing the cross-species integrated matrix. Print=TRUE prints out the cross-species integrated matrix and an RDS file of the Seurat object that results from the Seurat integration workflow.
#'
#' @export
#'
#' @examples

cross_species_integration <- function(matrix_path, min_cells = 3, min_features = 200, k_weight=100, n_pcs=30, UMAP_dim = c(1:30), FindNeighbors_dim=c(1:30), FindClusters_res = 0.5,  print=TRUE){
  
  
  ### create Seurat object
  
  #read in combined species matrix aqhumor_combined_mtx
  mat <- Read10X(matrix_path)
  #qc: filter for min cells and min features
  seurat.obj <- CreateSeuratObject(counts = mat, project="cross_species_integration", min.cells = min_cells, min.features=min_features)
  
  #qc may have filtered out cells. Remove missing cells from mat
  #seurat cells - the cells in the seurat object that survived filtering min.cells=3, min.features=200
  seurat_cells <- seurat.obj@assays[["RNA"]] %>% colnames()
  filtered_mat <- mat[,which(colnames(mat) %in% seurat_cells)]
  
  
  ### Annotate Seurat object with species, sample, cluster id,  cell type.
  
  #parse out annotations from the cell ids
  species_id <- sapply(strsplit(colnames(filtered_mat), "_"), function(v){return(v[1])})
  sample_id <- sapply(strsplit(colnames(filtered_mat), "_"), function(v){return(v[2])})
  cluster_id <- sapply(strsplit(colnames(filtered_mat), "_"), function(v){return(v[3])})
  cell_type_id <- sapply(strsplit(colnames(filtered_mat), "_"), function(v){return(v[4])})
  broad_cell_type_id <- sapply(strsplit(cell_type_id, "-"), function(v){return(v[1])})
  
  #Add annotation to Seurat obj
  seurat.obj <- AddMetaData(seurat.obj, species_id, col.name="species_id")
  seurat.obj <- AddMetaData(seurat.obj, sample_id, col.name = "sample_id")
  seurat.obj <- AddMetaData(seurat.obj, cluster_id, col.name = "cluster_id")
  seurat.obj <- AddMetaData(seurat.obj, cell_type_id, col.name = "cell_type_id")
  seurat.obj <- AddMetaData(seurat.obj, broad_cell_type_id, col.name = "broad_cell_type_id")
  
  seurat.obj@meta.data
  
  
  ### Integrate
  
  seurat.obj.list <- SplitObject(seurat.obj, split.by = "species_id")
  
  
  features <- SelectIntegrationFeatures(object.list = seurat.obj.list)
  
  anchors <- FindIntegrationAnchors(object.list = seurat.obj.list, anchor.features = features)
  
  #repeatable no matter what object you choose.
  seurat.obj.example <- seurat.obj.list[[1]]
  all_genes <- seurat.obj.example@assays[["RNA"]] %>% row.names()
  
  combined <- IntegrateData(anchorset = anchors, k.weight=k_weight, features.to.integrate = all_genes)
  
  
  combined<- ScaleData(combined, verbose=FALSE)
  combined<- RunPCA(combined, npcs=n_pcs, verbose=FALSE)
  combined<- RunUMAP(combined, reduction = "pca", dims=UMAP_dim)
  combined<- FindNeighbors(combined, reduction = "pca", dims = FindNeighbors_dim)
  combined<- FindClusters(combined, resolution=FindClusters_res)
  
  integrated<- combined[["integrated"]]@scale.data
  integrated <- as.sparse(integrated)
  
  #Save data
  if (print==TRUE){
    #create a directory for the matrix if it doesn't already exist
    if(!dir.exists("matrix")){
      dir.create("matrix")
    }
    #create a directory for the cross-species integrated matrix
    if(!dir.exists("matrix/cross-species_integration")){
      dir.create("matrix/cross-species_integration")
    }
    saveRDS(combined, "cross_species_integrated.rds")
    write10xCounts("aqhumor_cross_species_integrated_mtx", integrated, version="3")
    #move the matrix into the folder
    file.rename("aqhumor_cross_species_integrated_mtx", "matrix/cross-species_integration/aqhumor_cross_species_integrated_mtx")
    
    
  }
  
  return(combined)
  
}

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

#' id_errored_runs
#'
#' When `contml` runs error, it produces an empty `outtree` and `outfile`. We remove these, but now the index between the `contml` output files and scjack matrices do not match, leading to downstream analysis issues.
#' This function makes a list of matrices that produced an error when `contml` was run. We will then use this list to remove these matrices and proceed with downstream analyses.
#' This function is part of the workflow of the scjackknife analysis.
#'
#'
#' @importFrom utils write.table
#' 
#' @param outtrees_dir Path to directory with contml `outtrees`. This directory contains all `outtree` files created by contml, including for runs that errored (ie unprocessed outtrees directory)
#' @param matrix_name_pattern Prefix for matrix name. This function assumes that your matrix names end in `index_num_mtx` eg. for `means_5_reps_540_cells_20PC_after_PCA_index_mtx`, the matrix_name_pattern is `means_5_reps_540_cells_20PC_after_PCA`
#'
#' @return prints out a list of scjackknife matrix names that produced empty tree files when run through contml.
#' @export
#'
#' @examples
#'
id_errored_runs <- function(outtrees_dir, matrix_name_pattern){
  
  files_list <- list.files(path = outtrees_dir, full.names=TRUE, recursive=FALSE )
  
  empty_paths <- files_list[file.size(files_list) == 0]
  empty_files <- sapply(strsplit(empty_paths, "/"), function(v){return(v[length(v)])})
  
  #convert to matrices dir names
  names <- sapply(strsplit(empty_files, "_"), function(v){return(v[length(v)-1])})
  matrix_names <- paste0(matrix_name_pattern,"_", names,"_mtx")
  
  write.table(matrix_names, file = "errored_runs.txt", quote=FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)
}

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
#' @param file_name_prefix The prefix for the outtree file names
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
#' @param file_name File name for the labelled tree this function produces.
#' @param print Print out the tree file of the labelled tree. Tree is given in newick format.
#'
#' @return A phylo object of the labelled tree.
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


#' make_infile
#'
#' Converts a matrix to contml (PHYLIP) infile format.
#'
#' @importFrom Seurat Read10X
#' @importFrom utils write.table

#' @param matrix Matrix to create contml infile from. Use this if you have the matrix as an R object.
#' @param matrix_path Path to matrix to create contml infile from. Use this if you have a printed matrix in 10x feature-barcode format.
#' @param print Print the infile
#'
#' @return An object representing a matrix formatted in contml infile format. Contml requires a printed infile so print=TRUE is recommended.
#' @export
#'
#' @examples
make_infile = function(matrix, matrix_path, print){
  
  if (!missing(matrix_path) && missing(matrix)){
    mat <- Read10X(matrix_path)
  }
  
  if (!missing(matrix) && missing(matrix_path)){
    mat <- matrix
  }
  
  ntaxa <- nrow(mat)
  nchar <- ncol(mat)
  NA_col <- rep(NA, times=nrow(mat))
  phylip=""
  phylip <- cbind(1:nrow(mat), NA_col, NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col,NA_col, mat)
  phylip_final <- rbind(c(rep(NA, times=4), ntaxa, rep(NA, times=3), nchar, rep(NA, times=ncol(phylip)-9)), phylip)
  
  if (print == TRUE){
    print(write.table(phylip_final, file="infile", quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE, na= ""))}
  
  return(phylip_final)
  
}

#' make_pc_sweep_infiles
#'
#' This function produces the input files for the PC sweep experiment in Mah & Dunn (2023) (Fig.2).
#' It takes in the 919 cell 919 PC matrix (Fig. S1, step 7) and subsets the matrix to an increasing number of PCs. These matrix subsets are printed out as `contml` infiles.
#' If `calculate_PCA = TRUE` the function will first perform PCA on the matrix before subsetting and creating the infiles.
#'
#' @param matrix_path Path to a matrix that has had PCA performed on it.
#' @param run_PCA Boolean. If TRUE, calculate PCA on the input matrix before subsetting and producing infiles (i.e. the input matrix has not had PCA performed on it yet)
#' @param PC_range A range indicating the number of PCs to subset from the input matrix.
#'
#' @return This function prints out the `contml` infiles for the PC sweep experiment
#' @export
#'
#' @examples
#'
make_pc_sweep_infiles <- function(matrix_path, run_PCA = FALSE, PC_range=c(3:100)){
  
  
  mat <- Read10X(matrix_path) %>% as.matrix()
  
  if(run_PCA == TRUE){
    mat <- t(mat)
    
    #pca
    pc <- prcomp(mat, center=TRUE, scale=FALSE)
    x <- pc$x
    #normalize sd
    sds <- apply(x,2,sd)
    x.scaled <- t(t(x)/sds)
    apply(x.scaled, 2, sd)
    mat <- x.scaled
  }
  #set up directories to deposit infiles into
  if(!dir.exists("PC_sweep")){
    dir.create("PC_sweep")
  }
  
  #create a directory for the combined matrix
  if(!dir.exists("PC_sweep/infiles")){
    dir.create("PC_sweep/infiles")
  }
  
  
  lapply(PC_range, function(c){
    x.sweep.matrix <- mat[,1:c]
    #format to phylip
    phylip <- make_infile(matrix = x.sweep.matrix, print=FALSE)
    #write table
    write.table(phylip, file=paste0("PC_sweep/infiles/infile_pc_sweep_",c), quote=FALSE, sep=" ",   row.names=FALSE, col.names=FALSE, na= "")
    print(paste0("infile ", c))
  } #close function
  ) #close lapply
  
}

#' make_scjack_sample
#'
#' Create resampled matrices that contain 1 replicate cell per cell type group from a larger matrix that contains multiple replicates. This is used for running the scjackknife analysis.
#' Previously named make_scjack_sample54
#'
#' @importFrom utils write.table
#' @importFrom DropletUtils write10xCounts
#' @importFrom Seurat Read10X
#'
#' @param n_samples Number of scjackknife samples to create.
#' @param subsample_matrix_path Path to the matrix to sample from.
#' @param full_matrix_path Path to full unsubsetted matrix. PCA should already have been performed and PCs subsetted.
#' @return pca_matrix List of resampled single jackknife matrices.
#' @export
#'
#' @examples
make_scjack_sample <- function(n_samples, subsample_matrix_path, full_matrix_path, file_name){
  #load the full 919 cell 20 PC matrix for subsetting down below
  full_pca_mat <- Read10X(full_matrix_path)
  
  #make directories to organize
  if(!dir.exists("scjackknife")){
    dir.create("scjackknife")
  }
  
  if(!dir.exists(paste0("scjackknife/input_", file_name))){
    dir.create(paste0("scjackknife/input_", file_name))
  }
  
  if(!dir.exists(paste0("scjackknife/input_", file_name,"/scjack_selected_cells"))){
    dir.create(paste0("scjackknife/input_", file_name,"/scjack_selected_cells"))
  }
  
  if(!dir.exists(paste0("scjackknife/input_", file_name,"/scjack_matrices"))){
    dir.create(paste0("scjackknife/input_", file_name,"/scjack_matrices"))
  }
  
  if(!dir.exists(paste0("scjackknife/input_", file_name, "/scjack_infiles"))){
    dir.create(paste0("scjackknife/input_", file_name, "/scjack_infiles"))
  }
  
  sapply(n_samples, function(n_samples){
    
    #draw random subset of cells. 1 cell per cell type cluster per speciesu
    sampled_cells <- subset_combined_matrix(matrix_path = subsample_matrix_path, sample_size=1, print=FALSE)
    #record list of selected cells
    subset_list <- colnames(sampled_cells)
    
    #rank reduce matrix
    #pca_matrix <- run_pca(matrix_path = matrix_path, n_PCs = 20, subset = subset_list, print=FALSE)
    pca_matrix <- full_pca_mat[rownames(full_pca_mat) %in% subset_list,]
    
    #make infile
    phylip <- make_infile(matrix = pca_matrix, print = FALSE)
    
    ## print out scjackknife input files
    ## print outside of the two above functions so can have sequential file names
    
    ## record the subsampled cell ids
    write.table(subset_list, paste0("scjackknife_", file_name,"_selected_cells_",n_samples, ".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    file.rename(paste0("scjackknife_", file_name,"_selected_cells_",n_samples, ".txt"), paste0("scjackknife/input_", file_name,"/scjack_selected_cells/scjackknife_", file_name,"_selected_cells_",n_samples, ".txt"))
    ## write out subasmpled matrix
    write10xCounts(paste0("scjackknife_", file_name,"_",n_samples, "_mtx"),  pca_matrix, version="3")
    file.rename(paste0("scjackknife_", file_name,"_",n_samples, "_mtx"), paste0("scjackknife/input_", file_name,"/scjack_matrices/scjackknife_", file_name,"_",n_samples, "_mtx"))
    ## phylip infile
    write.table(phylip, file=paste0("infile_scjackknife_", file_name,"_", n_samples), quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE, na= "")
    file.rename(paste0("infile_scjackknife_", file_name,"_", n_samples), paste0("scjackknife/input_", file_name,"/scjack_infiles/infile_scjackknife_", file_name,"_", n_samples))
    
    return(pca_matrix)
    
  } #close function
  
  )#close sapply
  
} #close scjackknife_make_input()

#' plot_eigs
#'
#' Perform PCA and produce eigenvalue plots for a matrix. These plots are 1) eigenvalue vs PC, 2) eigenvalue cumulative sum, 3) eigenvalue vs percent variance, and 4) percent variance vs PC.
#'
#' @import ggplot2
#' @importFrom Seurat Read10X
#' @importFrom stats prcomp sd
#' @importFrom factoextra get_eigenvalue
#' @importFrom grDevices pdf dev.off
#'
#' @param matrix_path Path to matrix PCA is to be performed on.
#' @param highlight The number of the PC you would like to highlihgt in the plot. The dot representing this PC will be highlighted in red. Default is PC20, used in Mah & Dunn (2023)
#' @param print Boolean. Print out plots.
#' @return plots A list containing the four eigenvalue plots as ggplot objects
#' @export
#'
#' @examples
plot_eigs <- function(matrix_path, highlight=20, print=TRUE){
  
  #within-species normalized, combined and integrated cross-species. NOT subsetted.
  #unfiltered raw data
  mat <- Read10X(matrix_path)
  
  
  #transpose the matrix for the pca function.
  mat <- as.matrix(mat) %>% t()
  
  #pca
  pca <- prcomp(mat, center=TRUE, scale=FALSE)
  
  #eigenvalues, % var, cumulative % var
  eigvals <- get_eigenvalue(pca)
  
  #build plots
  #1. eigenvalue vs PC
  plot.df <- data.frame(order = 1:length(eigvals$eigenvalue), eigenvalue = eigvals$eigenvalue)
  plot.df.highlight <- plot.df[highlight,]
  p1 <- ggplot(plot.df, aes(x=order, y=eigenvalue)) + geom_point() + geom_point(data = plot.df.highlight, colour = "red" ) + theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  #2. eigenvalue cumulative sum
  plot.df <- data.frame(order = 1:length(cumsum(eigvals$eigenvalue)), eigenvalue_cum_sum = cumsum(eigvals$eigenvalue))
  p2 <- ggplot(plot.df, aes(x=order, y=eigenvalue_cum_sum)) + geom_point() +theme_bw() +xlab("order") +ylab("eigenvalue cumulative sum")
  #3 eigenvalue vs percent variance
  plot.df <- data.frame(eigenvalue = eigvals$eigenvalue, percent_variance = eigvals$variance.percent)
  p3 <- ggplot(plot.df, aes(x=eigenvalue, y=percent_variance)) + geom_point() +theme_bw() + xlab("eigenvalue") + ylab("percent variance")
  
  #how much variance is described by the first few PCs?
  #aqhumor_cross_species_integrated_mtx/: 26.65411
  print(paste0("The amount of variance described by the first ",highlight," PCs is: ",sum(plot.df$percent_variance[1:highlight])))
  
  # percent variance vs PC
  plot.df <- data.frame(component = 1:length(eigvals$variance.percent), percent_variance = eigvals$variance.percent)
  plot.df.highlight <- plot.df[highlight,]
  p4 <- ggplot(plot.df, aes(x = component, y = percent_variance)) + geom_point(size=0.5)  + xlab("Principal Component") + ylab("Percent Variance") + geom_point(data = plot.df.highlight, colour = "red", size=0.5) + theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_text(size=6.5), axis.text.y = element_text(size=6.5), axis.title = element_text(size = 7))
  
  if (print==TRUE){
    
    if(!dir.exists("eigenvalue_plots")){
      dir.create("eigenvalue_plots")
    }
    
    
    
    #4 variance vs
    #line plot of eigenvalues ordered
    pdf(file="eigenvalue_plots/1_eigenvalue_vs_PC.pdf", width=8, height=5)
    print(p1)
    dev.off()
    
    #cumulative sum
    pdf(file="eigenvalue_plots/2_eigenvalue_cumulative_sum.pdf", width=8, height=5)
    print(p2)
    dev.off()
    
    #variance of each component vs eigenvalue
    pdf(file="eigenvalue_plots/3_eigenvalue_vs_variance.pdf", width=5.5, height=5)
    print(p3)
    dev.off()
    
    #component variance
    pdf(file="eigenvalue_plots/4_eigenvalue_component_variance.pdf",  width=8, height=5)
    print(p4)
    dev.off()
  }
  
  plots <- list("1_eigenvalue_vs_PC" = p1, "2_eigenvalue_cumulative_sum" = p2, "3_eigenvalue_vs_variance" = p3, "4_eigenvalue_component_variance" = p4)
  
  return(plots)
  
}

#' plot_pc_sweep
#'
#' Plots the variables calculated by calc_pc_sweep_var. Produces Figs 2A, 2B, 2C. Default values are what were used in Mah & Dunn (2023).
#'
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#'
#' @param rank The range of the ranks you want to plot. This sets the x-axis for hte plots.The range must match the number of rows in the data frame of variables calculted by `calc_pc_sweep_var`.
#' @param var_df The dataframe of variable produced by cellphylo::calc_pc_sweep_var
#'
#' @return List containing plots for tree length, average edge length and star-ness score
#' @export
#'
#' @examples
plot_pc_sweep <- function(rank=c(3:100), var_df, print=TRUE){
  #assemble plot dfs
  plot.df <- data.frame(rank = rank, int.branch.len  = unlist(var_df$ave_int_edge), tip.branch.len = unlist(var_df$ave_tip_edge), star.measure = unlist(var_df$star_measure), TL = unlist(var_df$TL),sum.int.edge = unlist(var_df$sum_int_edge), sum.tip.edge = unlist(var_df$sum_tip_edge))
  
  #plot both on same plot - average branch length
  plot.df.2 <- data.frame(branch_lengths = c(plot.df$int.branch.len,plot.df$tip.branch.len), Edge_Type = c(rep("interior", length(plot.df$int.branch.len)), rep("tip", length(plot.df$tip.branch.len))), x = rep(rank, 2), check.names=FALSE)
  
  #plot sum of branch lengths instead of average
  plot.df.3 <- data.frame(branch_lengths = c(plot.df$sum.int.edge,plot.df$sum.tip.edge), Edge_Type = c(rep("interior", length(plot.df$int.branch.len)), rep("tip", length(plot.df$tip.branch.len))), x = rep(rank, 2), check.names=FALSE)
  
  
  #plots
  
  #TL
  tree_len.plot <- ggplot(data = plot.df, aes(x=rank, y=TL)) + geom_path() + geom_point() + xlab("PC") + ylab("Tree Length") + theme_bw() + scale_colour_grey(start = 0, end = .9)
  
  
  #ave edgelength
  #ave_edge_len.plot <- ggplot(data = plot.df.2, aes(x=x, y=branch_lengths, group=Edge_Type, colour=Edge_Type)) + geom_path() + geom_point(aes(shape=Edge_Type)) + xlab("PC") + ylab("Average Edge Length")  + theme_bw()  + theme(legend.position=c(0.80,0.5),  legend.box.background = element_rect(colour = "black"),legend.key.size = unit(0.25, "cm")) + scale_colour_grey(start = 0, end = .7)
  
  #sum of edge lengths
  edge_len.plot <- ggplot(data = plot.df.3, aes(x=x, y=branch_lengths, group=Edge_Type, colour=Edge_Type)) + geom_path() + geom_point(aes(shape=Edge_Type)) + xlab("PC") + ylab("Sum of Edge Lengths")  + theme_bw()  + theme(legend.position=c(0.80,0.5),  legend.box.background = element_rect(colour = "black"),legend.key.size = unit(0.25, "cm")) + scale_colour_grey(start = 0, end = .7)
  
  
  
  #star measure
  star.plot <- ggplot(data = plot.df, aes(x=rank, y=star.measure)) +  geom_path( ) + geom_point() + xlab("PC") + ylab("Star-ness") + theme_bw() + scale_colour_grey(start = 0, end = .9)
  
  if (print==TRUE){
    
    #create a PC sweep directory if doesn't already exist
    if(!dir.exists("PC_sweep")){
      dir.create("PC_sweep")
    }
    
    #create a directory for the combined matrix
    if(!dir.exists("PC_sweep/plots")){
      dir.create("PC_sweep/plots")
    }
    
    #print
    #TL
    pdf(file="1_tree_length.pdf", width=8, height=5)
    print(tree_len.plot)
    dev.off()
    
    #ave edge length
    #pdf(file="2_average_edge_length.pdf", width=8, height=5)
    #print(ave_edge_len.plot)
    #dev.off()
    
    #sum of edge lengths
    pdf(file="2_sum_edge_length.pdf", width=8, height=5)
    print(edge_len.plot)
    dev.off()
    
    #starness score
    pdf(file="3_star_measure.pdf", width=8, height=5)
    print(star.plot)
    dev.off()
    
    file.rename("1_tree_length.pdf", "PC_sweep/plots/1_tree_length.pdf")
    file.rename("2_sum_edge_length.pdf", "PC_sweep/plots/2_sum_edge_length.pdf")
    file.rename("3_star_measure.pdf", "PC_sweep/plots/3_star_measure.pdf")
    
    
  }
  
  return(list(tree_len.plot, edge_len.plot, star.plot))
  
}

#' run_pca
#'
#' Centres the value of a matrix around the mean, performs PCA on this matrix and normalizes variance by the standard deviation.
#' The number of principal components to retain in the final PCA matrix may be indicated with the "n_PCs" option.
#' It is recommended that PCA be performed on the full, integrated unsubsampled matrix. The final PCA matrix may subset to a list of cells specified by the user by the `"subset" or "subset_file_path" option.
#'
#' @importFrom Seurat Read10X as.sparse
#' @importFrom stats prcomp sd
#' @importFrom utils read.table
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to matrix to perform PCA on
#' @param n_PCs Number of principal components to reduce rank to
#' @param subset A list of cell ids to subsample from the input matrix and retain in the resulting PCA matrix this function creates. If no subset list is provided all cells are kept.
#' @param subset_file_path A text file containing a single column of cell ids to subsample from the input matrix and retain in the resulting PCA matrix this function creates.
#' @param print Print out the PCA matrix
#'
#' @return mat.sparse A matrix PCA has been performed on in sparse matrix format. This matrix may have subsetted from the full PCA matrix calculated from the input matrix by number of PCs and cell subset to retain.
#' @export
#'
#' @examples
run_pca <- function(matrix_path, n_PCs=20, subset, subset_file_path, print=TRUE){
  
  #read in matrix
  mat <- Read10X(matrix_path)
  
  #orient matrix
  mat <- as.matrix(mat) %>% t()
  
  #pca
  pc <- prcomp(mat,center=TRUE, scale=FALSE)
  x <- pc$x
  #normalize by sd after rotation
  sds <- apply(x,2,sd)
  x_norm <- t(t(x)/sds)
  apply(x_norm, 2, sd)
  
  #take just the first 20 pcs
  #check n PCs < n cells and > 1 PC. Matrix was rotated so n cells = n rows
  if (n_PCs > nrow(mat)){
    stop("The number of PCs to keep must be fewer than the number of cells")
  }
  if (n_PCs < 2){
    stop("The number of PCs to keep must be greater than 1.")
  }
  x_norm <- x_norm[,1:n_PCs]
  print(paste0("Subsetted to ", ncol(x_norm), " PCs."))
  
  # subset x to 1 cell per cell type cluster per species
  #Was a path or object given?
  if (!missing(subset_file_path) && missing(subset)){
    selected <- read.table(subset_file_path) %>% unlist()
  }
  
  if (!missing(subset) && missing(subset_file_path)){
    selected <- subset %>% unlist()
  }
  
  if (missing(subset) && missing(subset_file_path)){
    selected <- rownames(x_norm)
  }
  
  x.sub <- x_norm[rownames(x_norm) %in% selected,]
  print(paste0("Subsetted ", nrow(mat), " cells to ", nrow(x.sub), " cells."))
  
  #write matrix
  mat.sparse <- as.sparse(x.sub)
  
  if (print==TRUE){
    #if(!dir.exists("matrix")){
    #  dir.create("matrix")
    #}
    
    #create a directory for the combined matrix
    #if(!dir.exists("matrix/PCA")){
    #  dir.create("matrix/PCA")
    #}
    
    #write matrix
    write10xCounts(path=paste0("contml_", nrow(x.sub), "cell_subset_pca_var_norm_", n_PCs, "PC_mtx"), x= mat.sparse, version="3")
    #file.rename(paste0("contml_", nrow(x.sub), "cell_subset_pca_var_norm+", n_PCs, "PC_mtx"), paste0("matrix/PCA/contml_",nrow(x.sub), "cell_subset_pca_var_norm_", n_PCs, "PC_mtx"))
    
  } #close if print
  
  return(mat.sparse)
  
  
} #close func

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

#' subset_cells
#'
#' This function randomly samples n replicates per cell type group from a larger counts matrix. The input matrix must be in cellphylo format.
#'
#' @importFrom Seurat Read10X as.sparse
#' @importFrom dplyr filter
#' @importFrom utils write.table
#' @importFrom DropletUtils write10xCounts
#'
#' @param path_to_matrix Path to matrix you want to subsample from. Matrix must be in cellphylo format.
#' @param n_reps Number of replicates per cell type group to be in the subsampled matrix.
#'
#' @return mat.sub A sparse matrix created by random subsampling of the input matrix. If print = TRUE, prints out the cell ids of the sampled replicates (subset_cells_list) and the subsampeld matrix in 10x feature-barcode format.
#'
#' @examples
#'
#'
subset_cells <- function(path_to_matrix, n_reps, print=TRUE){
  
  print("Reading in matrix...")
  mat <- Read10X(path_to_matrix)
  print("Matrix loaded.")
  
  #species identity
  species <- sapply(strsplit(colnames(mat)[1], "_"), function(v){return(v[1])})
  #original full cell id
  orig_ids <- colnames(mat)
  #cluster ids (cluster num + cell type)
  cluster_ids <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[3])})
  #organize ids into a dataframe
  ids.df <- data.frame(orig_ids = orig_ids, cluster_id = cluster_ids)
  
  #how many cluster ids
  cluster_ids.unique <- unique(cluster_ids)
  
  subset_dfs <- lapply(cluster_ids.unique, function(x){
    
    ids.df.subset <- ids.df %>% filter(cluster_id == x)
    random_subset <- sample(ids.df.subset$orig_ids, size = n_reps, replace=FALSE)
    random_subset <- as.data.frame(random_subset)
    return(random_subset)
  })
  
  combined_subset_dfs <- do.call(rbind, subset_dfs)
  
  
  #Find the cell ids of the randomly selected cell barcodes in the orignal integrated matrix. Old approach used indices. Save ids.
  selected <- combined_subset_dfs$random_subset
  
  #subset matrix
  indices <- which(colnames(mat) %in% selected)
  mat.sub <-mat[,indices]
  
  #write matrix
  mat.sparse <- as.sparse(mat.sub)
  if(print==TRUE){
    
    if(!dir.exists("matrix")){
      dir.create("matrix")
    }
    
    #create a directory for the within-species normalized and integrated matrices
    if(!dir.exists("matrix/subset")){
      dir.create("matrix/subset")
    }
    
    if(!dir.exists("matrix/subset/subset_cells_list")){
      dir.create("matrix/subset/subset_cells_list")
    }
    
    write.table(selected, paste0("subset_cells_list_", species, ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
    file.rename( paste0("subset_cells_list_", species, ".txt"),  paste0("matrix/subset/subset_cells_list/subset_cells_list_", species, ".txt"))
    
    write10xCounts(paste0("matrix_subset_",n_reps,"_rep_", species), mat.sub, version="3")
    file.rename(paste0("matrix_subset_",n_reps,"_rep_", species),  paste0("matrix/subset/matrix_subset_",n_reps,"_rep_", species))
    
  }
  
  return(mat.sub)
  
} #close function

#' subsample_combined_matrix
#'
#' Randomly subsample a combined matrix containing cells from multiple species. Cell ids must be in cellphylo format.
#'
#' @importFrom Seurat Read10X as.sparse
#' @importFrom dplyr filter
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to a combined matrix of cells from multiple species.
#' @param sample_size The number of replicates per cell type group to be sampled
#' @param print Print out a file listing the cell ids of selected cells and the subsampled matrix in 10x format.
#'
#' @return mat.sparse The subsampled matrix in sparse matrix format.
#' @export
#'
#' @examples
#'
subset_combined_matrix <- function(matrix_path, sample_size, print){
  
  #load matrix
  mat <- Read10X(matrix_path)
  #parse cell ids
  orig_ids <- colnames(mat)
  cluster_ids <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[3])})
  species_id <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[1])})
  species_cluster_ids <- paste0(species_id, "_", cluster_ids)
  
  #organize annotations into a dataframe
  ids.df <-data.frame(orig_ids = orig_ids, cluster_ids = cluster_ids, species_cluster_ids = species_cluster_ids)
  species_cluster_ids.unique <- unique(species_cluster_ids)
  length(species_cluster_ids.unique)
  
  
  subset_dfs <- lapply(species_cluster_ids.unique, function(x){
    
    ids.df.subset <- ids.df %>% filter(species_cluster_ids== x)
    random_subset <- sample(ids.df.subset$orig_ids, size = sample_size, replace=FALSE)
    random_subset <- as.data.frame(random_subset)
    return(random_subset)
  })
  
  combined_subset_dfs <- do.call(rbind, subset_dfs)
  n_subset <- length(species_cluster_ids.unique)*sample_size
  
  
  #Find the indices of the randomly selected cell barcodes in the orignal integrated matrix
  selected <- combined_subset_dfs$random_subset
  indices <- which(colnames(mat) %in% selected)
  
  #subset matrix
  mat.sub <-mat[,indices]
  mat.sparse <- as.sparse(mat.sub)
  
  if (print==TRUE){
    
    #write list of subsampled cells
    write.table(selected, paste0("aqhumor_cross_species_",n_subset, "_subset.txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    
    #write subsampled matrix
    write10xCounts(paste0("aqhumor_cross_species_",n_subset, "_subset_mtx"), mat.sparse, version="3")
  }
  
  return(mat.sparse)
  
}

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' within_species_integration
#'
#' This function normalizes counts using a negative binomial model (SCTransform) and performs within-species batch integration using the Seurat workflow. Matrices must be in cellphyo format.
#'
#' @import Seurat
#' @importFrom DropletUtils write10xCounts
#'
#' @param path_to_matrix Path to a species' matrix. Cell ids must be in cellphylo format.
#' @param n_features The number of features to use for integration using Seurat::SelectIntegrationFeatures. Default is 3000, used in Mah & Dunn 2023.
#' @param k_weight The k.weight to use with Seurat::IntegrateData. This must be less than or equal to the number of cells in the smallest batch. Default is 100, used in Mah & Dunn 2023.
#' @param print Boolean. Print out the Seurat object as an RDS file and extract and print out the pearson residuals, corrected UMI counts and integrated pearson residuals from the Seurat object. These files will be deposited in `matrix/within-species_analysis`.
#'
#' @return combined.sct A Seurat object featuring the normalized and batch-integrated matrices.
#' If print=TRUE, it will also print out:
#' An RDS file of the Seurat object and 3 matrices in 10x format featuring values from the normalized and batch-integrated Seurat object: 1. sct_pearson_residuals (SCT assay "scale.data" slot), 2; sct_corrected_UMI_counts (SCT assay "counts" slot), 3.sct_integrated_pearson_residuals (integrated assay "scale.data" slot)
#'
#' @export
#'
#' @examples
#'

within_species_integration <- function(path_to_matrix, n_features=3000, k_weight = 100, print=TRUE){
  
  
  #load counts matrix
  mat <- Read10X(path_to_matrix)
  
  #min.features=200 from basic clustering tutorial
  seurat.obj <- CreateSeuratObject(counts = mat, project = "within_species", min.features=200)
  
  #parse annotations from cell names
  cell_type_id <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[4])})
  sample_id <-  sapply(strsplit(colnames(mat), "_"), function(v){return(v[2])})
  species <- sapply(strsplit(colnames(mat), "_"), function(v){return(v[1])})
  #species name for making files
  species_name <- species[1]
  
  #annotate seurat object
  seurat.obj <- AddMetaData(seurat.obj, sample_id, col.name="sample")   #sample aka "batch"
  seurat.obj <- AddMetaData(seurat.obj, cell_type_id , col.name="cell_type")
  seurat.obj <- AddMetaData(seurat.obj, species , col.name="species")
  
  #split by sample type
  seurat.obj.list <- SplitObject(seurat.obj, split.by="sample")
  
  #normalize with SCTransform
  #seurat.obj.list <- lapply(X = seurat.obj.list, FUN = function(x){ x<- SCTransform(x, return.only.var.genes = FALSE, residual.features = NULL, min_cells = 1, conserve.memory = FALSE)})
  seurat.obj.list <- lapply(X = seurat.obj.list, FUN = function(x){ x<- SCTransform(x)})
  
  features <- SelectIntegrationFeatures(object.list = seurat.obj.list, nfeatures=n_features)
  seurat.obj.list <- PrepSCTIntegration(object.list = seurat.obj.list, anchor.features = features)
  
  anchors <- FindIntegrationAnchors(object.list = seurat.obj.list, normalization.method = "SCT", anchor.features = features)
  
  #make character vector of all gene names
  seurat.obj.example <- seurat.obj.list[[1]]
  all_genes <- seurat.obj.example@assays[["RNA"]] %>% row.names()
  
  #if samples have less than 100 cells must set k.weight < 100
  combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight=k_weight, features.to.integrate = all_genes)
  
  #save matrices
  #SCTransform pearson residuals
  sct_pearson_residuals <- combined.sct[["SCT"]]@scale.data
  sct_pearson_residuals <- as.sparse(sct_pearson_residuals)
  #lognormalized, SCT-corrected UMI counts
  sct_corrected_UMI_counts <- combined.sct[["SCT"]]@counts
  sct_corrected_UMI_counts <- as.sparse(sct_corrected_UMI_counts )
  #integrated pearson residuals
  sct_integrated_pearson_residuals <- combined.sct[["integrated"]]@scale.data
  sct_integrated_pearson_residuals <- as.sparse(sct_integrated_pearson_residuals)
  
  if (print==TRUE){
    #save object
    saveRDS(combined.sct, paste0("within_species_integrated_",species_name,".rds"))
    
    #write matrices
    write10xCounts(paste0("matrix_sct_pearson_residuals_", species_name), sct_pearson_residuals, version="3")
    write10xCounts(paste0("matrix_sct_corrected_UMIs_", species_name), sct_corrected_UMI_counts, version="3")
    write10xCounts(paste0("matrix_sct_integrated_pearson_residuals_", species_name), sct_integrated_pearson_residuals, version="3")
    
    #create a matrix directory if doesn't already exist
    if(!dir.exists("matrix")){
      dir.create("matrix")
    }
    
    #create a directory for the within-species normalized and integrated matrices
    if(!dir.exists("matrix/within-species_norm_and_integration")){
      dir.create("matrix/within-species_norm_and_integration")
    }
    
    if(!dir.exists("matrix/within-species_norm_and_integration/pearson_residuals")){
      dir.create("matrix/within-species_norm_and_integration/pearson_residuals")
    }
    
    if(!dir.exists("matrix/within-species_norm_and_integration/corrected_UMIs")){
      dir.create("matrix/within-species_norm_and_integration/corrected_UMIs")
    }
    
    if(!dir.exists("matrix/within-species_norm_and_integration/integrated_pearson_residuals")){
      dir.create("matrix/within-species_norm_and_integration/integrated_pearson_residuals")
    }
    #move matrices into the directory
    file.rename(paste0("matrix_sct_pearson_residuals_", species_name), paste0("matrix/within-species_norm_and_integration/pearson_residuals/matrix_sct_pearson_residuals_",species_name))
    file.rename(paste0("matrix_sct_corrected_UMIs_", species_name), paste0("matrix/within-species_norm_and_integration/corrected_UMIs/matrix_sct_corrected_UMIs_", species_name))
    file.rename(paste0("matrix_sct_integrated_pearson_residuals_", species_name), paste0("matrix/within-species_norm_and_integration/integrated_pearson_residuals/matrix_sct_integrated_pearson_residuals_",species_name))
    
  }
  
  #return the normalized and integrated Seurat object
  return(combined.sct)
  
}

#' wrangle_matrix
#'
#' Converts cell identifiers from van Zyl et al. (2020)'s format to cellphylo format.
#'
#' @importFrom Seurat Read10X
#' @importFrom DropletUtils write10xCounts
#' @importFrom utils read.table
#' @importFrom dplyr mutate left_join
#'
#' @param path_to_mat Path to matrix to be wrangled into cellphylo format. Matrix must be in 10x format.
#' @param metadata_file Path to wrangled metadata file from van Zyl et al. (2020)
#' @param species Name of species the matrix belongs to
#'
#' @return mat A wrangled matrix with cell ids in cellphylo format
#' @export
#'
#' @examples
wrangle_matrix <- function(path_to_mat, metadata_file, species){
  
  #load wrangled meta data file
  meta_data <- read.table(file=metadata_file, header=TRUE, sep="\t")
  
  #load matrix
  mat <- Read10X(path_to_mat)
  #extract cell ids in order they are found in matrix
  matrix_cell_ids <- colnames(mat) %>% as.data.frame()
  colnames(matrix_cell_ids) <- "cell_identifier"
  
  #reorder meta data df to match mat colnames
  names.df <- left_join(matrix_cell_ids, meta_data, by= "cell_identifier")
  
  #parse out meta data info
  sample_names <- names.df$sample_id
  cell_barcode_ids <- names.df$cell_barcode
  cluster_id <- names.df$cluster_id
  cell_type_id <- gsub('[0-9]+', "",names.df$cluster_id)
  species_id <- species
  
  #create new cell names
  names.df <- names.df %>% mutate(new_cell_ids = paste0(species_id, "_", sample_names,"_",cluster_id,"_", cell_type_id,"_", cell_barcode_ids))
  
  #replace matrix names
  colnames(mat) <- names.df$new_cell_ids
  
  #write the new matrix
  write10xCounts(paste0("matrix_wrangled_",species), mat, version="3")
  
  #create a folder for matrices, if not already present
  if(!dir.exists("matrix")){
    dir.create("matrix")
  }
  
  if(!dir.exists("matrix/wrangled")){
    dir.create("matrix/wrangled")
  }
  
  file.rename(paste0("matrix_wrangled_",species), paste0("matrix/wrangled/matrix_wrangled_",species))
  
  return(mat)
  
}

#' wrangle_metadata
#'
#' This function takes in the meta data file `all_five_species_metafile.csv` provided by van Zyl et al. (2020) and parses out the cell identifier, cluster id, sample id and cell barcode.
#'
#' @importFrom dplyr mutate select
#' @importFrom utils read.csv
#' @param path_to_meta_file A path to the metadata file provided by van Zyl et al. (2020).
#'
#' @return meta_data A data frame with columns for cell identifier, cluster id, sample id and cell barcode parsed from the meta data file.
#' @export
#'
#' @examples
#'
#'
#'
wrangle_metadata <- function(path_to_meta_file){
  #library(tidyverse)
  
  #read in meta data file
  meta_data <- read.csv(path_to_meta_file, header=TRUE)
  #format - remove first row
  meta_data <- meta_data[-1,]
  
  #remove all `'` and spaces - messes up data frame actions
  meta_data$Cluster <- gsub("\\s", "_", meta_data$Cluster)
  meta_data$Cluster <- gsub("'", "", meta_data$Cluster)
  
  #save original cluster names
  meta_data <- meta_data %>% mutate(orig_cluster_label = meta_data$Cluster)
  colnames(meta_data) <- c("cell_identifier", "cluster", "orig_cluster_label")
  
  #cell barcodes in metafile end in -1 but in matrix end in .1
  #reserve "-" for cell type label delimeter
  meta_data$cell_identifier<- gsub("-1", ".1", meta_data$cell_identifier)
  
  #replace cell names
  meta_data$cluster <- gsub("NK/T","NKT", meta_data$cluster)
  meta_data$cluster  <- gsub("NK/T_cell","NKT", meta_data$cluster)
  meta_data$cluster  <- gsub("NKT_cell","NKT", meta_data$cluster)
  meta_data$cluster  <- gsub("Nonmyelinating_schwann_cell","SchwannCell-nmy", meta_data$cluster)
  meta_data$cluster <- gsub("Myelinating_schwann_cell","SchwannCell-my", meta_data$cluster)
  meta_data$cluster <- gsub("Schwann_cell","SchwannCell", meta_data$cluster)
  meta_data$cluster  <- gsub("BeamCella", "Beam-A", meta_data$cluster)
  meta_data$cluster  <- gsub("Beam_A", "Beam-A", meta_data$cluster)
  meta_data$cluster  <- gsub("BeamA", "Beam-A", meta_data$cluster)
  meta_data$cluster  <- gsub("BeamCellb","Beam-B", meta_data$cluster)
  meta_data$cluster  <- gsub("Beam_X","Beam-X", meta_data$cluster)
  meta_data$cluster <- gsub("Beam_Y","Beam-Y", meta_data$cluster)
  meta_data$cluster <- gsub("Ciliary_muscle", "CiliaryMuscle", meta_data$cluster)
  meta_data$cluster <- gsub("CribiformJCT","JCT", meta_data$cluster)
  meta_data$cluster <- gsub("B_cell","BCell", meta_data$cluster)
  meta_data$cluster <- gsub("VascularEndo","Endo-Vasc", meta_data$cluster)
  meta_data$cluster <- gsub("Vascular_Endothelium","Endo-Vasc", meta_data$cluster)
  meta_data$cluster <- gsub("ScEndo","Endo-Schlemms", meta_data$cluster)
  meta_data$cluster <- gsub("Corneal_endothelium","Endo-Corneal", meta_data$cluster)
  meta_data$cluster <- gsub("Endothelium","Endo", meta_data$cluster)
  meta_data$cluster <- gsub("CornealEpi","Epi-Corneal", meta_data$cluster)
  meta_data$cluster <- gsub("Corneal_epithelium","Epi-Corneal", meta_data$cluster)
  meta_data$cluster <- gsub("Pigmented_ciliary_epithelium","Epi-CiliaryPigment", meta_data$cluster)
  meta_data$cluster <- gsub("Nonpigmented_ciliary_epithelium","Epi-CiliaryNonPigment", meta_data$cluster)
  meta_data$cluster <- gsub("Pigmented_epithelium","Epi-Pigment", meta_data$cluster)
  meta_data$cluster <- gsub("SChlemms_Canal","SchlemmsCanal", meta_data$cluster) #typo fixed
  meta_data$cluster <- gsub("Collector_channel","CollectorChn", meta_data$cluster)
  
  #cluster id
  #collapse all underscores eg. 15_Pericyte -> 15Pericyte. Humans have no cluster number.
  cluster_id <- sub("_", "", meta_data$cluster)
  
  #sample ID
  sample_id <- sapply(strsplit(meta_data$cell_identifier, "_"), function(v){return(v[1])})
  cell_barcode <- sapply(strsplit(meta_data$cell_identifier, "_"), function(v){return(v[2])})
  
  #add annotations to meta data table
  meta_data <- meta_data  %>% mutate(cluster_id = cluster_id, sample_id = sample_id, cell_barcode = cell_barcode)
  meta_data <- meta_data %>% select(cell_identifier, cluster_id, sample_id, cell_barcode)
  #print out the wrangled meta data file
  utils::write.table(meta_data, "all_five_species_metafile_wrangled.csv", quote=FALSE, sep = "\t", col.names=TRUE, row.names=FALSE)
  
  if(!dir.exists("ann")){
    dir.create("ann")
  }
  
  file.rename("all_five_species_metafile_wrangled.csv", "ann/all_five_species_metafile_wrangled.csv")
  
  
  return(meta_data)
  
  
}

#' id_errored_runs_mean
#'
#' @param outtrees_dir ath to directory with contml `outtrees`. This directory contains all `outtree` files created by contml, including for runs that errored (ie unprocessed outtrees directory)
#' @param matrix_name_pattern Prefix for matrix name. This function assumes that your matrix names end in `index_num_mtx` eg. for `scjackknife_540_TBE_1_mtx`, the matrix_name_pattern is `scjackknife_540_TBE`
#'
#' @return NULL. Prints out a table with matrix ids that produced errored runs
#' @export
#'
#' @examples
id_errored_runs_mean <- function(outtrees_dir, matrix_name_pattern){
  
  files_list <- list.files(path = outtrees_dir, full.names=TRUE, recursive=FALSE )
  
  empty_paths <- files_list[file.size(files_list) == 0]
  empty_files <- sapply(strsplit(empty_paths, "/"), function(v){return(v[length(v)])})
  
  #convert to matrices dir names
  names <- sapply(strsplit(empty_files, "_"), function(v){return(v[length(v)-1])})
  ### CHANGED from cellphylo2 fun id-errored-runs
  #matrix_names <- paste0(matrix_name_pattern,"_", names,"_mtx")
  matrix_names <- paste0(matrix_name_pattern,"_", names)
  
  write.table(matrix_names, file = "errored_runs_mean.txt", quote=FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)
}

#' make_scjack_sample_mean
#' 
#' Create resampled matrices where each value is the average of 5 replicate cells per cell type group from a larger matrix. This is used for running the scjackknife analysis for Fig. 5 (averaged tree)
#' 
#' @param n_samples Number of scjackknife samples to create
#' @param matrix_path Path to matrix to sample from. Cell ids must be in cellphylo format
#' @param sample_size Number of cells to sample and then average per cell type group
#' @param limit_cells If you want to limit which cell ids to select, provide the list of cell ids to subsample from. Provide a file with a single column of cell ids.
#' @param file_name String to append to scjack file names
#' @param print Print out selected cell ids, averaged matrices and rownames of matrices, infiles for each scjackknife sample
#'
#' @return NULL. This function is meant to print out scjack samples and does not return an object.
#' @export
#'
#' @examples
make_scjack_sample_mean <- function(n_samples, matrix_path, sample_size, limit_cells, file_name, print){
  
  #make directories to organize
  if(!dir.exists("average")){
    dir.create("average")
  }
  
  if(!dir.exists("average/scjackknife")){
    dir.create("average/scjackknife")
  }
  
  if(!dir.exists(paste0("average/scjackknife/input_", file_name))){
    dir.create(paste0("average/scjackknife/input_", file_name))
  }
  
  if(!dir.exists(paste0("average/scjackknife/input_", file_name,"/scjack_selected_cells"))){
    dir.create(paste0("average/scjackknife/input_", file_name,"/scjack_selected_cells"))
  }
  
  if(!dir.exists(paste0("average/scjackknife/input_", file_name,"/scjack_matrices"))){
    dir.create(paste0("average/scjackknife/input_", file_name,"/scjack_matrices"))
  }
  
  if(!dir.exists(paste0("average/scjackknife/input_", file_name, "/scjack_infiles"))){
    dir.create(paste0("average/scjackknife/input_", file_name, "/scjack_infiles"))
  }
  
  if(!dir.exists(paste0("average/scjackknife/input_", file_name, "/scjack_matrix_cell_ids"))){
    dir.create(paste0("average/scjackknife/input_", file_name, "/scjack_matrix_cell_ids"))
  }
  
  sapply(1:n_samples, function(n){
    
    sub_means_mat <- subset_mean_matrix(matrix_path = matrix_path, sample_size=sample_size, print=print, file_name=file_name, limit_cells=limit_cells)
    phylip <- make_infile(matrix = sub_means_mat, print = print)
    
    #rename files and put into folders
    #matrix
    file.rename(paste0("means_",sample_size, "_reps_", file_name),paste0("average/scjackknife/input_", file_name,"/scjack_matrices/means_",sample_size, "_reps_", file_name, "_", n))
    #selected cells
    file.rename(paste0("selected_cells_",sample_size, "_reps_", file_name,".txt"),paste0("average/scjackknife/input_", file_name,"/scjack_selected_cells/selected_cells_",sample_size, "_reps_", file_name,"_",n,".txt"))
    #matrix cell ids 
    file.rename(paste0("rownames_means_",sample_size, "_reps_",file_name,".txt"),paste0("average/scjackknife/input_", file_name, "/scjack_matrix_cell_ids/rownames_means_",sample_size, "_reps_",file_name,"_",n,".txt"))
    #infiles
    file.rename("infile",paste0("average/scjackknife/input_", file_name, "/scjack_infiles/infile_",sample_size, "_reps_",file_name,"_",n))
    
  })
  
  
}

#' subset_mean_matrix
#' 
#' Randomly subsample a combined matrix and average each subsample for each cell type group. Cell ids must be in cellphylo format.
#'
#' @importFrom Seurat Read10X as.sparse
#' @importFrom DropletUtils write10xCounts
#' @importFrom dplyr filter
#' @importFrom utils write.table
#' 
#' @param matrix_path Path to matrix to sample from. Cell ids must be in cellphylo format
#' @param sample_size Number of cells to sample and then average per cell type group
#' @param print BOOLEAN. Print out selected cell ids, averaged matrices and rownames of matrices
#' @param file_name String to append to output matrix name.
#' @param limit_cells List of cell ids. If you want to limit which cell ids to select, provide the list of cell ids to subsample from.
#' @return means_mat.sparse The subsampled and averaged matrix in sparse matrix format.
#' @export
#'
#' @examples
subset_mean_matrix <- function(matrix_path, sample_size, print, file_name, limit_cells){
  #already made when making mean mat
  #919 cells - 10 cells per cell type cluster, 92 cell type clusters
  #matrix_path = "../matrix_final/integrated/cross_species/cross_species_integrated/aqhumor_cross_species_integrated_mtx/"
  #perform PCA
  #n_PCs = 20
  #run_pca(matrix_path, n_PCs, print=TRUE)
  
  #PARAMETER mat
  #load in new pca mat
  mat <- Read10X(matrix_path)
  mat <- t(as.matrix(mat))
  
  #PARAMETER limit_cells
  #if no list of cells to limit to, use all cell ids
  if (missing(limit_cells)){  
    limit_cells <- colnames(mat)
  }
  
  if (!missing(limit_cells)){  
    limit_cells <- unlist(limit_cells)
  }
  
  #parse out annotations
  orig_ids <- limit_cells
  species <- sapply(strsplit(orig_ids, "_"), function(v){return(v[1])})
  cluster_id <- sapply(strsplit(orig_ids, "_"), function(v){return(v[3])})
  cell_type_id <- sapply(strsplit(orig_ids, "_"), function(v){return(v[4])})
  species_cluster <- paste0(species, "_", cluster_id)
  #92 species cluster ids
  unique.species_cluster <- unique(species_cluster)
  
  #organize into id.df
  ids.df <-data.frame(orig_ids = orig_ids, cluster_ids = cluster_id, species_cluster_ids = species_cluster)
  
  #PARAMETER sample_size
  #sample_size=5
  
  #same as means but subsample to 5 cells then mean
  means <- lapply(1:length(unique.species_cluster), function(i){
    
    #subset ids.df
    ids.df.subset <- ids.df %>% filter(species_cluster_ids == unique.species_cluster[i])
    
    #randomly 5 cells
    random_subset <- sample(ids.df.subset$orig_ids, size = sample_size, replace=FALSE)
    #random_subset <- as.data.frame(random_subset)
    
    #subset mat
    mat.subset <- mat[,which(colnames(mat) %in% random_subset)]
    
    #keep track of which cells were selected
    selected_cells <-  colnames(mat.subset)
    
    #average the selected cells
    means.df <- rowMeans(mat.subset) %>% as.data.frame()
    colnames(means.df) <- unique.species_cluster[i]
    
    return(list(means.df = means.df, selected_cells = selected_cells))
    
  }) #close apply
  
  
  #PARAMETER print BOOLEAN
  if(print==TRUE)
  {
    #extract selected_cells from means list
    selected <- lapply(means,'[[',"selected_cells") %>% unlist()
    
    #print out selected cells
    write.table(selected, paste0("selected_cells_",sample_size, "_reps_", file_name,".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  means_mat <- lapply(means,'[[',"means.df")
  means_mat <- as.data.frame(means_mat)
  #dashes had been replaced with ".". Not sure why, make dashes again
  
  fixed_names  <- gsub("\\.", "-", colnames(means_mat))
  species <- sapply(strsplit(fixed_names, "_"), function(v){return(v[1])})
  cluster <-sapply(strsplit(fixed_names, "_"), function(v){return(v[2])})
  cell_type <- gsub('[0-9+]', "", cluster)
  new_names <- paste0(species, "_sample_", cluster, "_", cell_type, "_barcode")
  colnames(means_mat) <- new_names
  rownames(means_mat) <- rownames(mat)
  
  #return back to post PCA orientation
  means_mat <- t(means_mat)
  
  means_mat.sparse <- as.sparse(means_mat)
  
  if(print==TRUE){
    #print matrix
    write10xCounts(path=paste0("means_",sample_size, "_reps_", file_name), means_mat.sparse, version="3")
    
    #print rownames of matrix
    row_names <- rownames(means_mat.sparse)
    write.table(row_names, paste0("rownames_means_",sample_size, "_reps_",file_name,".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  return(means_mat.sparse)
  
}
