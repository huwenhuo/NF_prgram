#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters for Mouse EM-Seq ----
params.samplesheet    = "${projectDir}/samplesheet_mouse_emseq_s32_38.csv"
params.contrast_sheet = "${projectDir}/samplesheet_mouse_emseq_s32_38.csv"
params.output         = "results_emseq"
params.genome         = "mm10_sw"

// Import modules
include { TRIM_FASTQ } from './modules/trim_fastq'
include { BISMARK    } from './modules/align_dna'       // Or update module name accordingly
include { MULTIQC    } from './modules/multiqc'

workflow {
    ch_input_meta = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            def genome_info      = params.genomes[params.genome]
            
            // Using sampleID (or fallback to sample column) as unique gsm_id
            meta.gsm_id          = row.sampleID
            meta.sample_name     = row.sample_name
            meta.group           = row.genotype
            meta.bismark_index   = genome_info.bismark_index
            meta.fasta           = genome_info.fasta
            
            def r1 = file(row.fastq_1)
            def r2 = row.fastq_2 ? file(row.fastq_2) : null
            meta.r1 = r1
            meta.r2 = r2
            meta.mode = r2 ? "PE" : "SE"
            return meta
        }

    ch_trimmed = TRIM_FASTQ(ch_input_meta)

    ch_align_input = ch_trimmed.map { meta, trimmed_files ->
        def m = meta.clone()
        def trim_list = trimmed_files instanceof List ? trimmed_files : [trimmed_files]
        m.trim_r1 = trim_list[0]
        m.trim_r2 = trim_list.size() > 1 ? trim_list[1] : null
        return m
    }

    ch_aligned = BISMARK(ch_align_input)
}