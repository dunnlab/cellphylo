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
#' @print Boolean. Print out the Seurat object as an RDS file and extract and print out the pearson residuals, corrected UMI counts and integrated pearson residuals from the Seurat object. These files will be deposited in `matrix/within-species_analysis`.
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
