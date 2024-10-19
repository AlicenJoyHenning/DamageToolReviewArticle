#!/bin/bash

# CONTEXT 
# For running ensembleKQC, specific metrics are required to be inputted in a csv format.
# These are generated through the first script, extractFeatures.py. The outputs
# are then fed into the second script, runEnsembleKQC.py, to generate the damaged labels.

# Define directories
input_dir="/Users/alicen/Projects/ReviewArticle/input/"
intermediate_dir="/Users/alicen/Projects/ReviewArticle/EnsembleKQC_input/"
output_dir="/Users/alicen/Projects/ReviewArticle/EnsembleKQC_output/"

# Loop through all .csv files in the input directory
for input_file in "$input_dir"*.csv; do

  # Extract the base name of the file 
  base_name=$(basename "$input_file" | cut -d'_' -f1)
  
  # Define the intermediate file name
  intermediate_file="${intermediate_dir}${base_name}_input.csv"
  
  # Run the first Python script
  python ./extractFeatures.py "$input_file" "$intermediate_file" human --normalize=false
  
  # Define the output file name
  output_file="${output_dir}${base_name}_input.csv"
  
  # Run the second Python script
  python ./runEnsembleKQC.py --input_path="$intermediate_file" --labeled=false --output_path="$output_file"
done