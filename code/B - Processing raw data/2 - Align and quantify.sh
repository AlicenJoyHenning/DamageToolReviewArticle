#! /bin/bash

# Script to perform STARsolo alignment 
# Edit the < > enclosed fields to match your file paths and names. 
# For the FASTQ files, the following shows the case where  more than one lane is present, and the paths are separate with commas. If only one lane is present, use single input. 
# Note: This script is the general case (3' protocol). For alternatives, see the Walk through document

#PBS -P <Project>
#PBS -N Align_and_Quantify
#PBS -l select=1:ncpus=24:mem=100gb
#PBS -l walltime=04:00:00
#PBS -q serial
#PBS -M <email>

# Change to the directory where the script is located
cd $PBS_O_WORKDIR
eval "$(conda shell.bash hook)"

# Input to STARsolo function
genome=<star_index>
read2=<sample_L001_R2.fastq.gz,sample_L002_R2.fastq.gz,sample_L003_R2.fastq.gz,sample_L004_R2.fastq.gz>
read1=<sample_L001_R1.fastq.gz,sample_L002_R1.fastq.gz,sample_L003_R1.fastq.gz,sample_L004_R1.fastq.gz>
whitelist=<3M-february-2018.txt/737K-august-2016.txt>
output=<string_for_prefix_of_output_files> 

# Activate virtual environment
conda activate <STAR_env_path>

STAR --runThreadN 20 \
     --genomeDir "$genome" \
     --readFilesIn "$read2" "$read1" \
     --soloType Droplet \
     --soloCBwhitelist "$whitelist" \
     --readFilesCommand gunzip -c \
     --soloCBlen 16 \
     --soloUMIlen 12 \
     --outFileNamePrefix "$output" \
     --soloCellFilter EmptyDrops_CR \
     --soloFeatures Gene Velocyto \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes CB


conda deactivate 

