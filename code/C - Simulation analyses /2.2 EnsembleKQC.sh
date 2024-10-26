#!/bin/bash
 
# Like section B, this script runs EnsembleKQC using two scripts, extractFeatures.py, 
# followed by runEnsembleKQC.py, to generate the damaged labels for the simulated data.

# Define directories
input_dir="/Users/alicen/Projects/ReviewArticle/simulated_ddqc_input/"
intermediate_dir="/Users/alicen/Projects/ReviewArticle/EnsembleKQC_simulated_input/"
output_dir="/Users/alicen/Projects/ReviewArticle/EnsembleKQC_simulated_output/"

# Loop through all .csv files in the input directory
for input_file in "$input_dir"*.csv; do

  # Extract the base name of the file (everything before the 4th _)
  base_name=$(basename "$input_file" | cut -d'_' -f1-4)
  
  # Define the intermediate file name
  intermediate_file="${intermediate_dir}${base_name}_input.csv"
  
  # Run the first Python script
  python ./extractFeatures.py "$input_file" "$intermediate_file" human --normalize=false
  
  # Define the output file name
  output_file="${output_dir}${base_name}_input.csv"
  
  # Run the second Python script
  python ./runEnsembleKQC.py --input_path="$intermediate_file" --labeled=false --output_path="$output_file"
done