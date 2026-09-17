#!/bin/bash

module load nextflow/24.10.0
module load singularity/3.9.9

nextflow run /work/InternalMedicine/s184335/sc/repos/NF_prgram/sept16_hsa.nf \
    -resume \
    -c /work/InternalMedicine/s184335/sc/repos/NF_prgram/nextflow.config \
    -w /work/InternalMedicine/s184335/sc/workspace/sisi/admera_hsa_sep16 \
    --genome GRCh38_sw \
    --contrast_sheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/samplesheet_human_rnaseq_s21_31.csv \
    --samplesheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/samplesheet_human_rnaseq_s21_31.csv \
    --output /work/InternalMedicine/s184335/sc/workspace/sisi/admera_hsa_sep16/results
