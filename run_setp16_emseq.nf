#!/bin/bash

module load nextflow/24.10.0
module load singularity/3.9.9

nextflow run /work/InternalMedicine/s184335/sc/repos/NF_prgram/sept16_emseq.nf \
    -resume \
    -c /work/InternalMedicine/s184335/sc/repos/NF_prgram/nextflow.config \
    -w /work/InternalMedicine/s184335/sc/workspace/sisi/admera_emseq_sep16 \
    --genome mm10_sw \
    --samplesheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/samplesheet_mouse_emseq_s32_38.csv \
    --output /work/InternalMedicine/s184335/sc/workspace/sisi/admera_emseq_sep16/results