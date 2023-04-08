#' cross_species_integration
#'
#' @import Seurat
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to combined cross-species, unintegrated matrix produced by cellphylo::combine_matrices. Cell ids must be in cellphylo format.
#'
#' @return combined The Seurat object containing the cross-species integrated matrix.
#' @export
#'
#' @examples
cross_species_integration <- function(matrix_path){


  ### create Seurat object

  #read in combined species matrix aqhumor_combined_mtx
  mat <- Read10X(matrix_path)
  #qc: filter for min cells and min features
  seurat.obj <- CreateSeuratObject(counts = mat, project="cross_species_integration", min.cells = 3, min.features=200)

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

  seurat.obj.example <- seurat.obj.list[[1]]
  all_genes <- seurat.obj.example@assays[["RNA"]] %>% row.names()

  combined <- IntegrateData(anchorset = anchors, k.weight=100, features.to.integrate = all_genes)


  ### Visualize

  combined <- ScaleData(combined, verbose=FALSE)
  combined <- RunPCA(combined, npcs=30, verbose=FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims=1:30)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
  combined <- FindClusters(combined, resolution=0.5)

  #Save data
  saveRDS(combined, paste0("cross_species_integrated.rds"))
  integrated<- combined[["integrated"]]@scale.data
  integrated <- as.sparse(integrated)
  write10xCounts("aqhumor_cross_species_integrated_mtx", integrated, version="3")

  return(combined)

}
