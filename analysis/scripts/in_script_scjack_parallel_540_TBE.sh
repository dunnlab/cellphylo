#!/bin/bash

#in_script_scjack_parallel_540_TBE.sh

#run contml on infile_n
#must not have an 'infile' indirectory - will replace with infile_n
#must have an (empty) outfile and outtree so can set

n=ITERATION_INDEX

cd /contml

echo "Run: ${n}"

#create random seed for jumble
RAND_NUM=$(shuf -i 1-2000000000 -n 1)
random_seed=$(($RAND_NUM * 2 + 1))
echo "seed: ${random_seed}"

outfile_new="outfile_scjack_540_TBE_${n}_${random_seed}"
outtree_new="outtree_scjack_540_TBE_${n}_${random_seed}"

infile="./scjack_infiles/infile_scjackknife_540_TBE_${n}"

echo "outfile: ${outfile_new}"
echo "outtree: ${outtree_new}"
echo "infile: ${infile}"

printf "${infile}\nF\n${outfile_new}\nC\nJ\n${random_seed}\n100\nG\nY\nF\n${outtree_new}\n"| ./phylip/phylip-3.697/exe/contml
#test
#printf "${infile}\nF\n${outfile_new}\nC\nJ\n${random_seed}\n2\nG\nY\nF\n${outtree_new}\n"| ./phylip/phylip-3.697/exe/contml

mv ${outtree_new} ./outtrees_scjack_540_TBE/${outtree_new}.tre
mv ${outfile_new} ./outfiles_scjack_540_TBE/${outfile_new}.txt


echo "Run ${n} complete"