# cellphylo

This repo contains code meant to facilitate reproduction of the analyses in Mah & Dunn (2023).  

# Replicate analyses
First run the file containing all functions used in the analyses: `cellphylo_functions.R`. SessionInfo can be found in `sessionInfo_R.txt`.  

All analyses can be replicated in full by running these notesbooks in the `analysis` directory, in order:  

>1_Wrangle_Data.Rmd: Format matrices into `cellphylo` format  
>2_Create_Matrices.Rmd: Create the within and cross-species integrated matrices (Fig. S1 steps 1-6)  
>3_PC_sweep_analysis.Rmd: Perform the PC sweep analysis (Fig. 2)  
>4_Jumble_analysis.Rmd: Calculate jumble scores (Fig. 4)  
>5_scJackknife_analysis.Rmd: Calculate scjackknife TBE scores (Fig. 4)  

# Directory Map  
Several files have been provided in the `analysis` directory to help with reproduction. These include certain input files, some helpful intermediate files and some output files. A thorough description of each of these files is provided in the notebooks as they walk through the analysis.  

A brief description of folders of `analysis`:  

analysis/  
 	ann/	Annotation files  
	jumble/  Output files for the jumble analysis  
	matrix/  
		cross-species_integration/	Combined, integrated matrix (Fig. S1, step 6)   
		subset/	Subsetted matrives for each species (Fig. S1, step 4)  
	PC_sweep/	Input and output files for the PC sweep analysis (Fig. 2)  
	scjackknife/	Input and output files for calculating scJackknife scores  
		booster/	Input files for `booster`, which was used to calculate TBE scores (Lemoine et al. 2018)  
		input/	Input files used to create the scjackknife trees  
		output/	Output files produced during the scjackknife analysis  
	scripts/	Bash scripts used to run `contml` for each analysis    
		jumble/	Scripts for jumble analysis  
		scjackknife/	Scripts for scJackknife analysis  
	tree_inference/	Matrices that were used to infer trees featured in the paper  
		54_cell_tree/	Matrix used to infer the 54 cell tree in Fig. 4  
		92_cell_tree/	Matrix used to infer the 92 cell tree in Fig. S4  


# Citation
Lemoine, F. et al. Renewing Felsenstein’s phylogenetic bootstrap in the era of big data. Nature 556, 452–456 (2018).
