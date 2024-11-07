#! /bin/bash

# Convert BAM to FASTQ files with Cellranger 
# In some cases, data is only available in BAM format. To convert to FASTQ, as required for alignment, 
# download the CellRanger softwar : https://www.10xgenomics.com/support/software/cell-ranger/downloads 
# The path to the downloaded (uncompressed) cellranger-v#-#-# must be specified inside the <download_path> below.

#PBS -P HEAL1360
#PBS -N BAM_FASTQ
#PBS -l select=1:ncpus=24:mem=100gb
#PBS -l walltime=04:00:00
#PBS -q serial
#PBS -M 27741931@sun.ac.za

# Change to the directory where the script is located
cd $PBS_O_WORKDIR

# Inputs 
bam_files=<path_to_BAM>

# Conversion
<download_path>/cellranger/cellranger-8.0.1/lib/bin/bamtofastq --nthreads=8 "bam_files" 

