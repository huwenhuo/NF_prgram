#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "${projectDir}/gse150984_samplesheet.tsv"
params.output      = "results"
params.genome      = "GRCh38_sw"

// Import processes from modules.nf in current working directory
include { 
    GET_SRR_ID;
    DOWNLOAD_FASTQ;
    TRIM_FASTQ; 
    STAR_TEALIGNMENT; 
    SC_TE;
    SC_TELOCAL;
    IRFINDER_FASTQ;
    MULTIQC
} from './modules.nf'

workflow {

    // 1. Build initial meta map channel from TSV samplesheet
    ch_input_meta = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [:]
            def genome_info    = params.genomes ? params.genomes[params.genome] : null
            meta.gsm_id        = row.gsm_id
            meta.sample_name   = row.sample
            meta.genome        = params.genome
            meta.star_index    = genome_info ? genome_info.star_index : null
            meta.bowtie2_index = genome_info ? genome_info.bowtie2_index : null
            return meta
        }

    // 2. Retrieve SRR IDs from GSM IDs
    ch_srr_out = GET_SRR_ID(ch_input_meta)

    // Print GSM ID and retrieved SRR ID directly to terminal
    ch_srr_out.view { meta, srr -> 
        "FOUND: GSM = ${meta.gsm_id} -> SRR = ${srr.trim()}" 
    }

    // Extract SRR ID and attach to meta map
    ch_srr_meta = ch_srr_out.map { meta, srr ->
        def updated_meta = meta.clone()
        def clean_srr = srr.toString().trim()
        updated_meta.srr_id = clean_srr
        return [updated_meta, clean_srr]
    }

    // 3. Download FASTQ files using retrieved SRR IDs
    ch_downloaded = DOWNLOAD_FASTQ(ch_srr_meta)
        .map { meta, files ->
            def file_list = files instanceof List ? files : [files]
            def updated_meta = meta.clone()

            updated_meta.r1   = file_list.find { it.name.contains('_1') } ?: file_list[0]
            updated_meta.r2   = file_list.find { it.name.contains('_2') }
            updated_meta.r3   = file_list.find { it.name.contains('_3') }
            updated_meta.mode = updated_meta.r2 ? (updated_meta.r3 ? "PE_plus_R3" : "PE") : "SE"

            return updated_meta
        }

    // 4. Trimming
    ch_trimmed = TRIM_FASTQ(ch_downloaded)

    // 5. Update meta with trimmed file locations
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

    // 6. Alignments
    ch_star_te_out = STAR_TEALIGNMENT(ch_aligned_input)

    // 7. Downstream TE analyses
    ch_te_input = ch_star_te_out.map { meta, bam, bai ->
        def m = meta.clone()
        def selected_genome = params.genomes[meta.genome]
        m.bam          = bam
        m.gtf          = selected_genome.gtf
        m.te_gtf       = selected_genome.te_gtf
        m.te_loc       = selected_genome.te_loc
        m.scTE_idx     = selected_genome.scTE_idx
        m.scTE_tx_idx  = selected_genome.scTE_tx_idx
        return m
    }

    SC_TE(ch_te_input)
    SC_TELOCAL(ch_te_input)

    // 8. Build inputs for IRFinder using trimmed FASTQ paths
    ch_irfinder_input = ch_trimmed.map { meta, trimmed_files ->
        def m = meta.clone()
        def selected_genome = params.genomes[meta.genome]
        def trim_list = trimmed_files instanceof List ? trimmed_files : [trimmed_files]
        if (meta.mode.startsWith("PE")) {
            m.trim_r1 = trim_list[0]
            m.trim_r2 = trim_list[1]
        } else {
            m.trim_r1 = trim_list[0]
            m.trim_r2 = null
        }
        m.irfinder_index = selected_genome.irfinder_index
        return m
    }
    
    IRFINDER_FASTQ(ch_irfinder_input)

    // 9. Collect QC logs for MultiQC
    ch_trim_qc = ch_trimmed.map { meta, files -> 
        (files instanceof List ? files : [files]).collect { it.getParent() }
    }.flatten()

    ch_star_qc = ch_star_te_out.map { meta, bam, bai -> 
        bam.getParent() 
    }

    ch_multiqc_inputs = ch_trim_qc
        .mix(ch_star_qc)
        .unique()
        .collect()

    MULTIQC(ch_multiqc_inputs)
}