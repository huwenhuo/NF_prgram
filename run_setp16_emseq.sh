#!/bin/bash
#SBATCH --partition=super     
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=200GB
#SBATCH --time=15-00:00:00
#SBATCH --job-name=emseq
#SBATCH --output=emseq.log

cd /work/InternalMedicine/s184335/sc/

source .bashrc

source ~/miniconda3/bin/activate ~/miniconda3/envs/r43

module load nextflow/24.10.0
module load singularity/3.9.9

nextflow run /work/InternalMedicine/s184335/sc/repos/NF_prgram/sept16_emseq.nf \
    -resume \
    -c /work/InternalMedicine/s184335/sc/repos/NF_prgram/nextflow.config \
    -w /work/InternalMedicine/s184335/sc/workspace/sisi/admera_emseq_sep16/work \
    --genome mm10_sw \
    --samplesheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/samplesheet_mouse_emseq_s32_38.csv \
    --output /work/InternalMedicine/s184335/sc/workspace/sisi/admera_emseq_sep16/results
