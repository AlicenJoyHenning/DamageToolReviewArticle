#! /bin/bash

# Building a STAR genome index
# Edit the < > enclosed fields to match your file paths and names

#PBS -P <Project>
#PBS -N Generate_Genome_Index
#PBS -l select=1:ncpus=24:mem=100gb
#PBS -l walltime=01:00:00
#PBS -q serial
#PBS -M <email>

# Change to the directory where the script is located
cd $PBS_O_WORKDIR
eval "$(conda shell.bash hook)"

# Input to STAR function 
output=<output_directory>
fasta=<genome.fa>
gtf=<genes.gtf>

# Load CHPC modules
module add chpc/python/anaconda/3-2021.11

# Activate STAR virtual environment 
conda activate <STAR_env_path>

# Run the index building 
STAR --runThreadN 23 \
     --runMode genomeGenerate \
     --genomeDir "$output" \
     --genomeFastaFiles "$fasta" \
     --sjdbGTFfile "$gtf" \

# Deactivate STAR virtual environment
conda deactivate 
