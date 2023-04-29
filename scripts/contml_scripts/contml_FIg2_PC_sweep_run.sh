#!/bin/bash

#pc sweep
#create a directory called `contml` and deposit a directory containing the PC sweep infiles into /contml/

cd /contml

mkdir outtrees_PC_sweep_919_final
mkdir outfiles_PC_sweep_919_final

for n in {3..100}
do

cd pc_sweep_infiles_final
cp infile_pca_var_norm_pc_sweep_${n} ../

cd ../
mv infile_pca_var_norm_pc_sweep_${n}  infile

#use species 1 as outgroup (default)
printf "C\nY\n"| ./phylip/phylip-3.697/exe/contml

mv outtree outtrees_PC_sweep_919_final/outtree_${n}_PC_sweep_919_final.tre
mv outfile outfiles_PC_sweep_919_final/outfile_${n}_PC_sweep_919_final.txt

echo "Run ${n}"


done

echo "PC_sweep 3-100 complete"

#leave docker container
exit

#shut down instance
sudo shutdown -h now