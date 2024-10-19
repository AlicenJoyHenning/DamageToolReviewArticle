#! /bin/bash

# Convert SRA to FASTQ files with sra toolkit
# In some cases, data is only available in SRA format. To convert to FASTQ, as required for alignment, 
# download the SRA toolkit : https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit. 
# The path of the downloaded sra_toolkit must be specified inside the <download_path> below.

#PBS -P <Project>
#PBS -N SRA_FASTQ
#PBS -l select=1:ncpus=24:mem=100gb
#PBS -l walltime=04:00:00
#PBS -q serial
#PBS -M <email>

# Change to the directory where the script is located
cd $PBS_O_WORKDIR

# Inputs 
sra_files=<path_to_SRA>

# Conversion 
<download_path>/sra_toolkit/sratoolkit.3.1.1-ubuntu64/bin/fastq-dump --split-files --gzip "$sra_files"

