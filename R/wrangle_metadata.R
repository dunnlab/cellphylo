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
