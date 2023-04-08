#' within_species_integr
#'
#' @import Seurat
#' @importFrom DropletUtils write10xCounts
#'
#' @param path_to_matrix Path to a species' matrix. Cell ids must be in cellphylo format.
#'
#' @return NULL RDS file of integrated Seurat object and 3 matrices in 10x format featuring values from the integrated Seurat object SCT assay "scale.data" slot (sct_pearson_residuals), SCT assay "counts" slot (sct_corrected_UMI_counts), and integrated assay "scale.data" slot (sct_integrated_pearson_residuals)
#' @export
#'
#' @examples
within_species_integr <- function(path_to_matrix){


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

  features <- SelectIntegrationFeatures(object.list = seurat.obj.list, nfeatures=3000)
  seurat.obj.list <- PrepSCTIntegration(object.list = seurat.obj.list, anchor.features = features)

  anchors <- FindIntegrationAnchors(object.list = seurat.obj.list, normalization.method = "SCT", anchor.features = features)

  #make character vector of all gene names
  seurat.obj.example <- seurat.obj.list[[1]]
  all_genes <- seurat.obj.example@assays[["RNA"]] %>% row.names()

  #if samples have less than 100 cells must set k.weight < 100
  combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight=100, features.to.integrate = all_genes)

  #save object
  saveRDS(combined.sct, paste0("combined_sct_",species_name,".rds"))

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

  #write matrices
  write10xCounts(paste0("sct_pearson_residuals_mtx_", species_name), sct_pearson_residuals, version="3")
  write10xCounts(paste0("sct_corrected_UMIs_mtx_", species_name), sct_corrected_UMI_counts, version="3")
  write10xCounts(paste0("sct_integrated_pearson_residuals_mtx_", species_name), sct_integrated_pearson_residuals, version="3")

}
