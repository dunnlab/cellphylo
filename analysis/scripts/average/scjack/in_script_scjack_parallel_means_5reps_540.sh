#!/bin/bash

#in_script_scjack_parallel_means_5reps_530.sh

#run contml on infile_n
#must not have an 'infile' in directory - will replace with infile_n
#must have an (empty) outfile and outtree so can set

n=ITERATION_INDEX

cd /contml

echo "Run: ${n}"

#create random seed for jumble
RAND_NUM=$(shuf -i 1-2000000000 -n 1)
random_seed=$(($RAND_NUM * 2 + 1))
echo "seed: ${random_seed}"

### CHANGE HERE
outfile_new="outfile_scjack_means_5reps_540_${n}_${random_seed}"
outtree_new="outtree_scjack_means_5reps_540_${n}_${random_seed}"

### CHANGE HERE
infile="./scjack_infiles/infile_5_reps_540_cells_20PC_after_PCA_${n}"

echo "outfile: ${outfile_new}"
echo "outtree: ${outtree_new}"
echo "infile: ${infile}"

#53 cell type groups: use G
#printf "${infile}\nF\n${outfile_new}\nC\nJ\n${random_seed}\n100\nG\nY\nF\n${outtree_new}\n"| ./phylip/phylip-3.697/exe/contml
#54 cell type groups: do not use G
printf "${infile}\nF\n${outfile_new}\nC\nJ\n${random_seed}\n100\nY\nF\n${outtree_new}\n"| ./phylip/phylip-3.697/exe/contml


### CHANGE HERE
mv ${outtree_new} ./outtrees_scjack_means_5reps_540/${outtree_new}.tre
mv ${outfile_new} ./outfiles_scjack_means_5reps_540/${outfile_new}.txt


echo "Run ${n} complete"