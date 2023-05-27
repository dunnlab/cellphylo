#!/bin/bash

#make_reps_scjack_means_5reps_530.s

### Modified from script produced by ChatGPT Dec 9 2022.
### Prompt: "Write a bash script that iteratively reads in a text file, copies and renames it a new iterative file name, and replaces a string in the file with the iteration index. Make it do it 100 times."

# Set the starting iteration index
index=1

### CHANGE HERE
# Set the filename of the original file
filename="in_script_scjack_parallel_means_5reps_540.sh"

# Set the string to be replaced in the file
replace="ITERATION_INDEX"

### CHANGE HERE
# Set the prefix for the new iterative filenames
new_filename_prefix="script_scjack_means_5reps_540_"

# Set the suffix for the new iterative filenames
new_filename_suffix=".sh"

# Set the maximum number of iterations
max_iterations=500


# Iterate until the original file does not exist or the maximum number of iterations is reached
while [ -f $filename ] && [ $index -le $max_iterations ]
do
  # Copy and rename the file with the current iteration index
  #don't need to do this
  #cp $filename $new_filename_prefix$index$new_filename_suffix

  # Replace the string in the copied file with the iteration index
  #sed -i "s/$replace/$index/g" $new_filename_prefix$index$new_filename_suffix
  #sed command did not work. Remove -i flag, input old file, pipe sed output to new file.
  sed "s/$replace/$index/g" $filename > $new_filename_prefix$index$new_filename_suffix
  
  #set permissions
  chmod u+x $new_filename_prefix$index$new_filename_suffix
  
  # Increment the iteration index
  index=$((index+1))
  
done
