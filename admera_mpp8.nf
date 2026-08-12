#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet.csv"
params.output      = "results"
params.genome      = "mm10" // Adjust default genome key as needed

selected_genome = params.genomes ? params.genomes[params.genome] : null

// Import processes from modules.nf
include { 
    TRIM_FASTQ; 
    ALIGN_RNA_STAR; 
    STAR_TEALIGNMENT; 
    TECOUNT; 
    TELOCAL; 
    IRFINDER_FASTQ;
    SC_TE;
    MULTIQC
} from './modules.nf'

workflow {

    // 1. Build initial meta map channel
    ch_input_meta = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            meta.gsm_id     = row.sample
            meta.r1         = row.fastq_1
            meta.r2         = row.fastq_2
            meta.mode       = row.fastq_2 ? "PE" : "SE"
            meta.star_index = selected_genome ? selected_genome.star_index : null
            return meta
        }

    // 2. Trimming
    ch_trimmed = TRIM_FASTQ(ch_input_meta)

    // 3. Update meta with trimmed file locations
    ch_aligned_input = ch_trimmed.map { meta, trimmed_files ->
        def updated_meta = meta.clone()
        if (meta.mode == "PE") {
            updated_meta.trim_r1 = trimmed_files[0]
            updated_meta.trim_r2 = trimmed_files[1]
        } else {
            updated_meta.trim_r1 = trimmed_files[0]
            updated_meta.trim_r2 = null
        }
        return updated_meta
    }

    // 4. Alignments
    ch_star_te_out = STAR_TEALIGNMENT(ch_aligned_input)

    // 5. Downstream TE analyses
    ch_te_input = ch_star_te_out.map { meta, bam, bai ->
        def m = meta.clone()
        m.bam      = bam
        m.gtf      = selected_genome.gtf
        m.te_gtf   = selected_genome.te_gtf
        m.te_loc   = selected_genome.te_loc
        m.scTE_idx = selected_genome.scTE_idx
        return m
    }

    TECOUNT(ch_te_input)
    TELOCAL(ch_te_input)
    SC_TE(ch_te_input)

    // 6. Build inputs for IRFinder using trimmed FASTQ paths and reference parameters
    ch_irfinder_input = ch_trimmed.map { meta, trimmed_files ->
        def m = meta.clone()
        if (meta.mode == "PE") {
            m.trim_r1 = trimmed_files[0]
            m.trim_r2 = trimmed_files[1]
        } else {
            m.trim_r1 = trimmed_files[0]
            m.trim_r2 = null
        }
        m.irfinder_index = selected_genome.irfinder_index
        m.img_irfinder   = params.img_irfinder
        return m
    }
    
    IRFINDER_FASTQ(ch_irfinder_input)

    // 7. Collect QC logs for MultiQC
    // Collect trimming files (e.g., fastp JSON/HTML) and alignment BAM outputs
    ch_trim_qc = ch_trimmed.map { meta, files -> files }.flatten()
    ch_star_qc = ch_star_te_out.map { meta, bam, bai -> bam }.flatten()

    // Combine all QC log paths into a single collection channel
    ch_multiqc_inputs = ch_trim_qc
        .mix(ch_star_qc)
        .collect()

    MULTIQC(ch_multiqc_inputs)
}