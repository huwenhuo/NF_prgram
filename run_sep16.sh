module load nextflow/24.10.0
module load singularity/3.9.9

nextflow run /work/InternalMedicine/s184335/sc/repos/NF_prgram/lsk_old.nf \
    -resume \
    -c /work/InternalMedicine/s184335/sc/repos/NF_prgram/nextflow.config \
    -w /work/InternalMedicine/s184335/sc/workspace/sisi/admera_sept16 \
    --genome mm10_sw \
    --contrast_sheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/samplesheet_mouse_rnaseq_s1_20.csv \
    --samplesheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/samplesheet_mouse_rnaseq_s1_20.csv \
    --output /work/InternalMedicine/s184335/sc/workspace/sisi/admera_sep16/results
