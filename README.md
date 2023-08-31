# cellphylo

This repo contains code meant to facilitate reproduction of the analyses in Mah & Dunn (2023), *Reconstructing cell type evolution across species through cell phylogenies of single-cell RNAseq data*.   

Link to bioRxiv:  https://www.biorxiv.org/content/10.1101/2023.05.18.541372v2   

# Replicate analyses
First run the file containing all functions used in the analyses: `cellphylo_functions.R`. SessionInfo can be found in `sessionInfo_R.txt`.  

All analyses can be replicated in full by running these notesbooks in the `analysis` directory, in order:  

*1_Wrangle_Data.Rmd*: Format matrices into `cellphylo` format  
*2_Create_Matrices.Rmd*: Create the within and cross-species integrated matrices (Fig. S1 steps 1-6)  
*3_PC_sweep_analysis.Rmd*: Perform the PC sweep analysis (Figure 2)
*4_Jumble_analysis.Rmd*: Infer the focal tree for Figure 4 and calculate jumble scores.  
*5_scJackknife_analysis.Rmd*: Calculate scjackknife scores for the tree in Figure 4.  
*6_Average_cells.Rmd*: Infer the focal tree for Figure 5 and calculate jumble and scjackknife scores.  
*7_GO_and_phenetic_trees.Rmd*: Perform the GO enrichment analysis (Table 1, S1, S6) and calculate the phenetic trees (Figs. S7-S10).  

# Directory Map  
Select input, intermediate and output files sufficient for reproduction of the analyses have been provided have in the `analysis` directory. A thorough description of each of these files is provided in the notebooks as they walk through the analysis.  

A brief description of folders of `cellphylo/analysis`:  

```   
analysis/
├── 1_Wrangle_Data.Rmd
├── 2_Create_Matrices.Rmd
├── 3_PC_sweep_analysis.Rmd
├── 4_Jumble_analysis.Rmd
├── 5_scJackknife_analysis.Rmd
├── 6_Average_cells.Rmd
├── 7_GO_and_phenetic_trees.Rmd
├── GO_enrichment
│   └── PCA_919x919_matrix.rds
├── PC_sweep
│   ├── input
│   │   └── pc_sweep_infiles_final.tar.gz
│   └── output
│       └── pc_sweep_919_outtrees.tar.gz
├── ann
│   ├── AqHumor_ensembl_orthologs.txt
│   └── all_five_species_metafile.csv
├── average
│   ├── 540_cell_subset_cell_ids.txt
│   ├── booster
│   │   ├── multiphylo_scjack_means_5reps_540_common_label.tre
│   │   ├── outtree_54mean_20PC_afterPCA_jumble_search_2525674179_common_label.tre
│   │   └── outtree_54mean_20PC_afterPCA_jumble_search_2525674179_common_label_best_tree_TBE.tre
│   ├── jumble
│   │   ├── input
│   │   │   ├── contml_54means_20PC_afterPCA_mtx
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── infile_54means_20PC_afterPCA
│   │   └── output
│   │       ├── 54mean_20PC_afterPCA_jumble_search_2525674179_scored_label.tre
│   │       └── outtrees_54mean_20PC_afterPCA_jumble_search_processed.tar.gz
│   └── scjackknife
│       ├── input
│       │   └── input_540_cells_20PC_after_PCA.tar.gz
│       └── output
│           └── output_processed
│               ├── errored_runs_mean.txt
│               ├── outtree_scjack_means_5reps_540_common_label.tar.gz
│               ├── outtrees_scjack_means_5reps_540_processed.tar.gz
│               └── scjack_matrices_processed.tar.gz
├── goseq_positive_BP_PC20_enriched.txt
├── jumble
│   ├── input
│   │   ├── contml_54cell_subset_pca_var_norm_20PC_mtx
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   └── infile_54_20PC
│   └── output
│       ├── outtree_54_jumble_search_1249177875_scored_label.tre
│       └── outtrees_54_jumble_search_processed.tar.gz
├── matrix
│   ├── cross-species_integration
│   │   ├── aqhumor_cross_species_92_subset.txt
│   │   └── aqhumor_cross_species_integrated_mtx
│   │       ├── barcodes.tsv.gz
│   │       ├── features.tsv.gz
│   │       └── matrix.mtx.gz
│   └── subset
│       ├── integrated_subset_mtx_human
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   └── matrix.mtx.gz
│       ├── integrated_subset_mtx_macF
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   └── matrix.mtx.gz
│       ├── integrated_subset_mtx_macM
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   └── matrix.mtx.gz
│       ├── integrated_subset_mtx_mouse
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   └── matrix.mtx.gz
│       └── integrated_subset_mtx_pig
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── phenetic_trees
│   ├── 54_means_pca_var_norm_20_PC_nj.tre
│   ├── 54_means_pca_var_norm_20_PC_upgma.tre
│   ├── 54_means_pca_var_norm_20_PC_wpgma_mcquitty.tre
│   └── NJ_scjack_analysis
│       ├── booster
│       │   └── multiphylo_nj_means_outtrees_scjack_540_TBE_common_label.tre
│       ├── output
│       │   ├── nj_means_outtrees_scjack_540_TBE_common_label.tar.gz
│       │   └── nj_trees_means_scjack_540.tar.gz
│       └── scored_trees
│           └── 54_means_BM_with_NJ_multiphylo.tre
├── scjackknife
│   ├── booster
│   │   ├── multiphylo_scjack_540_TBE_common_label.tre
│   │   ├── outtree_54_jumble_search_1249177875_common_label.tre
│   │   └── outtree_54_jumble_search_1249177875_scored_with_scjack540TBE.tre
│   ├── contml_919cell_subset_20PC_mtx
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── input
│   │   └── input_scjack_540_TBE.tar.gz
│   ├── matrix_540cell_cross_species_integrated_mtx
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── output
│       └── output_processed
│           ├── errored_runs.txt
│           ├── outtrees_scjack_540_TBE_common_label.tar.gz
│           ├── outtrees_scjack_540_TBE_processed.tar.gz
│           └── scjack_matrices_processed.tar.gz
├── scripts
│   ├── PC_sweep
│   │   └── contml_pc_sweep_919.sh
│   ├── average
│   │   ├── jumble
│   │   │   ├── run_jumble_parallel_54mean_20PC_afterPCA.sh
│   │   │   └── script_jumble_parallel_54mean_20PC_afterPCA.sh
│   │   ├── lsi
│   │   └── scjack
│   │       ├── in_script_scjack_parallel_means_5reps_540.sh
│   │       ├── make_reps_scjack_means_5reps_540.sh
│   │       ├── run_scjack_parallel_means_5reps_540.sh
│   │       └── scjack_scripts_mean.tar.gz
│   ├── jumble
│   │   ├── run_jumble_parallel_54.sh
│   │   └── script_jumble_parallel_54.sh
│   └── scjackknife
│       ├── in_script_scjack_parallel_540_TBE.sh
│       ├── make_reps_scjack_540_TBE.sh
│       ├── run_scjack_parallel_540_TBE.sh
│       └── scjack_scripts.tar.gz
└── tree_inference
    ├── 54_cell_tree
    │   └── contml_54cell_subset_pca_var_norm_20PC_mtx
    │       ├── barcodes.tsv.gz
    │       ├── features.tsv.gz
    │       └── matrix.mtx.gz
    └── 92_cell_tree
        ├── aqhumor_cross_species_92_subset.txt
        └── contml_92cell_subset_pca_var_norm_20PC_mtx
            ├── barcodes.tsv.gz
            ├── features.tsv.gz
            └── matrix.mtx.gz

54 directories, 94 files   

```    

# Supplementary Files  
Additional supplementary files referenced in the text are provided in `cellphylo/supplementary_files/` 

# Citations   
Mah, J.L., and Dunn, C.W. Reconstructing cell type evolution across species through cell phylogenies of single-cell RNAseq data. bioRxiv, (2023).  

Updated August 31 2023  

