module load nextflow/24.10.0

nextflow run /work/InternalMedicine/s184335/sc/repos/NF_prgram/gse150984.nf \
    -resume \
    -c /work/InternalMedicine/s184335/sc/repos/NF_prgram/nextflow.config \
    -w /work/InternalMedicine/s184335/sc/workspace/sisi/gse150984/work \
    --genome GRCh38_sw \
    --samplesheet /work/InternalMedicine/s184335/sc/repos/NF_prgram/gse150984_samplesheet.tsv \
    --output /work/InternalMedicine/s184335/sc/workspace/sisi/gse150984/work 

