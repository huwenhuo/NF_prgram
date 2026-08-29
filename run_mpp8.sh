module load nextflow/24.10.0
module load singularity/3.9.9

nextflow run /work/InternalMedicine/s184335/sc/repos/NF_prgram/mpp8.nf \
    -resume \
    -c /work/InternalMedicine/s184335/sc/repos/NF_prgram/nextflow.config \
    -w /work/InternalMedicine/s184335/sc/workspace/sisi/admera_mpp8/work \
    --genome mm10_sw \
    --contrast_sheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/mpp8_samplesheet_fixed.csv \
    --samplesheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/mpp8_samplesheet.csv \
    --output /work/InternalMedicine/s184335/sc/workspace/sisi/admera_mpp8/results

