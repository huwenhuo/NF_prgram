#!/bin/bash

TARGET_DIR="/archive/InternalMedicine/Chung_lab/shared/sc/sig_image"
export SINGULARITY_CACHEDIR="${TARGET_DIR}/.singularity_cache"
cd "$TARGET_DIR"

# Navigate to your target directory first
cd /archive/InternalMedicine/Chung_lab/shared/sc/sig_image

# Pull fastp 0.24.0
singularity pull fastp_0.24.0.sif docker://biocontainers/fastp:0.24.0--h8ee403a_0

# Pull MultiQC (latest)
singularity pull multiqc_latest.sif docker://ewels/multiqc:latest

# Pull SRA Tools (latest)
singularity pull -F /archive/InternalMedicine/Chung_lab/shared/sc/sig_image/sra-tools.sif docker://pegi3s/sratoolkit:3.1.0

# Pull Entrez direct tool
singularity pull entreztool.sif docker://pegi3s/entrez-direct:latest

# Pull Bowtie2 + SAMtools image
singularity pull bowtie2_2.4.1.sif docker://biocontainers/bowtie2:v2.4.1_cv1

# Pull samtools 1.24
singularity pull samtools_1.24.sif docker://staphb/samtools:1.24

# Pull STAR 2.7.11b from Docker Hub
singularity pull star_2.7.11b.sif docker://josousa/star:2.7.11b

# MACS3 3.0.3
singularity pull macs3_3.0.0b3.sif docker://joseespinosa/macs3:3.0.0b3

# Picard 3.4.0
singularity pull picard_3.4.0.sif https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0

# DeepTools 3.5.5
singularity pull deeptools_3.5.5.sif https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0

# Pull TEtranscripts
singularity pull tetranscripts.sif docker://mhammelllab/tetranscripts:latest

# Pull TElocal (latest)
singularity pull telocal_latest.sif docker://mhammelllab/telocal:latest

# Pull IRFinder
singularity pull irfinder_2.0.1.sif docker://cloxd/irfinder:2.0.1
