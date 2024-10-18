# From raw to processed count matrices

## 1. Create environment for future work 

STAR environment

```konsole
$ conda create -n STAR_environment star
```



## 2. Prepare input files for alignment

#### Helpful resource: https://sydney-informatics-hub.github.io/training-RNAseq/02-BuildAGenomeIndex/index.html

### 1.1 Reference genome Index
STAR aligns RNA-seq data to genomic DNA sequence. It does so in a splicing-aware manner, in that it accomodates for the natural “gaps” that occur when aligning RNA to genomic DNA sequence as a result of splicing. To run STAR, you need to have a **reference genome** to align to. For most species in the world, there is no reference and de novo assembly is required (hectic). But for many species, expecially model species like mice and macaques, references are readily available, _so_ readily available that you have a often have a choice...

For human and mouse data, the GENCODE annotations are the industry standard [source](https://sydney-informatics-hub.github.io/training-RNAseq/02-BuildAGenomeIndex/index.html). These are available on the [Gencode website](https://www.gencodegenes.org/). 

This code can be used for obtaining the genome sequence fasta files for the mouse. 
```konsole
  $ wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/GRCm39.primary_assembly.genome.fa.gz
```

In addition to the **reference genome**, you need the **genome annotations** that are associated with it. These can be obtained from the same location on the [Gencode website](https://www.gencodegenes.org/). For this, choose the option that says PRI ( contains the comprehensive gene annotation on the primary assembly (chromosomes and scaffolds) sequence regions).
```konsole
$ wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.primary_assembly.annotation.gtf.gz
```

Once the fasta and annotation files are downloaded, the genome index can be built. This can be done, among other ways, using a STAR command. But to use STAR the input .fasta and .gtf file need to be unzipped. This can be done in linux:
```konsole
$ gunzip file.fasta
$ gunzip file.gtf
```
Then, use the unzipped files as input in the following command where the run mode is specified to generate the genome (aka not the run mode for mapping reads). The output directory of the files is specified (genomeDir) as well as the paths for the fasta and gtf files.
```konsole
STAR --runThreadN 23 --runMode genomeGenerate --genomeDir /home/ahenning/lustre/work/Human/genome/index --genomeFastaFiles /home/ahenning/lustre/work/Human/genome/genome.fa --sjdbGTFfile /home/ahenning/lustre/work/Human/genome/annotation.gtf
```

 Note: this was run in an interactive session on the CHPC 
```konsole
$ qsub -I -l select=1:ncpus=24:mpiprocs=24:mem=100gb -P PROJECT_NAME -q serial -l walltime=2:00:00
$ conda activate STAR_environment
```


Alternatively, and preferably, you can generate the genome index in a job script and submit it to be run on the CHPC script (GenomeGenerate.sh) 

# 3. Alignment 

For scRNA-seq droplet based 10X protocols, there are 3 common cases that require different inputs and parameters during alignment. These are outlined here.

a) Chemistry v2 3' 

whitelist: 737K-august-2016.txt 

```konsole 

STAR --runThreadN 20 \
     --genomeDir "$genome" \
     --readFilesIn "$read2" "$read1 \
     --soloType Droplet \
     --soloCBwhitelist "$whitelist" \
     --readFilesCommand gunzip -c \
     --soloCBlen 16 \ 
     --soloUMIlen 10 \ 
     --outFileNamePrefix "$output" \
     --soloCellFilter EmptyDrops_CR \ 
     --soloFeatures Gene Velocyto \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes CM

```


b) Chemistry v2 5' "extra length"

whitelist: 737K-august-2016.txt 

STAR --runThreadN 20 \
     --genomeDir "$genome" \
     --readFilesIn "$read1" "$read2" \
     --soloType CM_UMI_Simple \
     --soloCBwhitelist "$whitelist" \
     --readFilesCommand gunzip -c \
     --soloBarcodeMate 1 \
     --clip5pNbases 39 0 \
     --soloCBstart 1 \
     --soloCBlen 16 \ 
     --soloUMIstart 17 \
     --soloUMIlen 10 \ 
     --outFileNamePrefix "$output" \
     --soloCellFilter EmptyDrops_CR \ 
     --soloFeatures Gene Velocyto \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes CM


c) Chemistry v3

whitelist: 3M-february-2018.txt

```konsole 

STAR --runThreadN 20 \
     --genomeDir "$genome" \
     --readFilesIn "$read2" "$read1 \
     --soloType Droplet \
     --soloCBwhitelist "$whitelist" \
     --readFilesCommand gunzip -c \
     --soloCBlen 16 \ 
     --soloUMIlen 12 \ 
     --outFileNamePrefix "$output" \
     --soloCellFilter EmptyDrops_CR \ 
     --soloFeatures Gene Velocyto \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes CM

```












