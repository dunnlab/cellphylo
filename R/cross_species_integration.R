#' cross_species_integration
#'
#' Integrate counts across species using the Seurat workflow. Default values are what were used in Mah & Dunn (2023)
#'
#' @import Seurat
#' @importFrom DropletUtils write10xCounts
#'
#' @param matrix_path Path to combined matrix. This combined matrix features cells from multiple species but is not yet cross-species integrated. This matrix was created using `cellphylo::combine_matrices`. Cell ids must be in cellphylo format.
#' @param min_cells  Specify the min number of cells a feature must be expressed in to be retained. This is the min.cells option of Seurat::CreateSeuratObject. Default is 3.
#' @param min_features Specaify the min number of features a cell must be expressed to be retained. This is the min.features option of Seurat::CreateSeuratObject. Default is 200.
#' @param k_weight The k.weight to use with Seurat::IntegrateData. This must be less than or equal to the number of cells in the smallest batch. Default is 100.
#' @param print Boolean. Print out the matrix and Seurat object of cross-species integrated data.
#'
#' @return combined The Seurat object containing the cross-species integrated matrix. Print=TRUE prints out the cross-species integrated matrix and an RDS file of the Seurat object that results from the Seurat integration workflow.
#'
#' @export
#'
#' @examples

cross_species_integration <- function(matrix_path, min_cells = 3, min_features=200, k_weight=100,  print=TRUE){


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

  seurat.obj.example <- seurat.obj.list[[1]]
  all_genes <- seurat.obj.example@assays[["RNA"]] %>% row.names()

  combined <- IntegrateData(anchorset = anchors, k.weight=k_weight, features.to.integrate = all_genes)


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
  saveRDS(combined, paste0("cross_species_integrated.rds"))
  integrated<- combined[["integrated"]]@scale.data
  integrated <- as.sparse(integrated)
  write10xCounts("matrix_cross-species_integrated", integrated, version="3")
  #move the matrix into the folder
  file.rename("matrix_cross-species_integrated", "matrix/cross-species_integration/matrix_cross-species_integrated")


  }

  return(combined)

}
