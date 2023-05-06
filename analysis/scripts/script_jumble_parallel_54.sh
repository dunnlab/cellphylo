#!/bin/bash

#Max Jumble seed is 10 digits (4294967293). Draw random number and make it odd x*2 +1. Max value of x is 2147483645.5.

cd /contml

RAND_NUM=$(shuf -i 1-2000000000 -n 1)
random_seed=$(($RAND_NUM * 2 + 1))
echo "Seed ${random_seed}"


#make file names
outfile_new="outfile_54_jumble_search_${random_seed}"
outtree_new="outtree_54_jumble_search_${random_seed}"
echo ${outfile_new}
echo ${outtree_new}

#select J with random seed, jumble 100x and select tree with highest likelihood
printf "F\n${outfile_new}\nC\nJ\n${random_seed}\n100\nG\nY\nF\n${outtree_new}\n"| ./phylip/phylip-3.697/exe/contml
#test
#printf "F\n${outfile_new}\nC\nJ\n${random_seed}\n2\nG\nY\nF\n${outtree_new}\n"| ./phylip/phylip-3.697/exe/contml

mv ${outtree_new} outtrees_54_jumble_search/${outtree_new}.tre
mv ${outfile_new} outfiles_54_jumble_search/${outfile_new}.txt


echo "jumble run complete"


