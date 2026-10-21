#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "${projectDir}/samplesheet_mouse_rnaseq_s1_20.csv"
params.output      = "results"
params.genome      = "mm10_sw"

// Import or define modules/processes
include { TRIM_FASTQ } from './modules/trim_fastq'
include { ALIGN_RNA_STAR } from './modules/align_rna_star.nf'

// ==========================================
// WORKFLOW: Trim & Align Only
// ==========================================
workflow {
    ch_input_meta = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            def genome_info = params.genomes ? params.genomes[params.genome] : null
            meta.gsm_id        = row.sampleID
            meta.sample_name   = row.sample_name
            meta.group1        = row.group
            meta.rep           = row.rep
            meta.genome        = params.genome
            meta.star_index    = genome_info ? genome_info.star_index : null
            
            def r1 = file(row.fastq_1)
            def r2 = row.fastq_2 ? file(row.fastq_2) : null
            meta.r1 = r1
            meta.r2 = r2
            meta.mode = r2 ? "PE" : "SE"

            return meta
        }

    // 1. Trim Fastq files
    ch_trimmed = TRIM_FASTQ(ch_input_meta)

    // 2. Format trimmed outputs back into meta for STAR alignment
    ch_aligned_input = ch_trimmed.map { meta, trimmed_files ->
        def updated_meta = meta.clone()
        def trim_list = trimmed_files instanceof List ? trimmed_files : [trimmed_files]
        if (meta.mode.startsWith("PE")) {
            updated_meta.trim_r1 = trim_list[0]
            updated_meta.trim_r2 = trim_list[1]
        } else {
            updated_meta.trim_r1 = trim_list[0]
            updated_meta.trim_r2 = null
        }
        return updated_meta
    }

    // 3. Run STAR Alignment with GeneCounts quantification
    ALIGN_RNA_STAR(ch_aligned_input)
}