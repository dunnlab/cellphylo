#!/bin/bash

#pc sweep
#create a directory called `contml` and deposit a directory `infiles/` containing the PC sweep infiles into /contml/

cd /contml

mkdir pc_sweep_outtrees
mkdir pc_sweep_outfiles

for n in {3..100}
do

cd infiles
cp infile_pc_sweep_${n} ../

cd ../
mv infile_pc_sweep_${n}  infile

#use species 1 as outgroup (default)
printf "C\nY\n"| ./phylip/phylip-3.697/exe/contml

mv outtree c_sweep_outtrees/outtree_pc_sweep_${n}tre
mv outfile pc_sweep_outfiles/outfile_pc_sweep_${n}.txt

echo "Run ${n}"


done

echo "PC_sweep 3-100 complete"

#leave docker container
exit

#shut down instance
sudo shutdown -h now