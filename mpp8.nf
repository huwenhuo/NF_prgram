#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "${projectDir}/mpp8_samplesheet.csv"
params.output      = "results"
params.genome      = "mm10_sw"

// Import processes from modules in current working directory
include { TRIM_FASTQ } from './modules/trim_fastq'
include { ALIGN_RNA_STAR } from './modules/align_rna_star'
include { STAR_TEALIGNMENT } from './modules/star_tealignment'
include { SAMTOOLS_INDEX } from './modules/samtools_index'

include { TECOUNT } from './modules/tecount'
include { MERGE_TECOUNTS } from './modules/merge_tecounts'
include { DESEQ2_TECOUNT } from './modules/deseq2_tecount'
include { VOLCANO_TECOUNT } from './modules/volcano_tecount'
include { HEATMAP_TECOUNT } from './modules/heatmap_tecount'

include { DESEQ2_CODING } from './modules/deseq2_coding'
include { HEATMAP_ANALYSIS } from './modules/heatmap_analysis'
include { PATHWAY_ANALYSIS } from './modules/pathway_analysis'
include { VOLCANO_PLOT } from './modules/volcano_plot'

include { TELOCAL } from './modules/telocal'
include { MERGE_TELOCAL } from './modules/merge_telocal'
include { EDGER_TELOCAL } from './modules/edger_telocal'

include { SC_TE } from './modules/sc_te'
include { SC_TELOCAL } from './modules/sc_telocal'
include { MERGE_SCTELOCAL } from './modules/merge_sctelocal'
// include { EDGER_SCTELOCAL } from './modules/edger_sctelocal'

include { IRFINDER_FASTQ } from './modules/irfinder_fastq'
include { MERGE_IRFINDER } from './modules/merge_irfinder'
include { DESEQ2_IRFINDER } from './modules/deseq2_irfinder'
include { HEATMAP_IRFINDER } from './modules/heatmap_irfinder'
include { VOLCANO_IRFINDER } from './modules/volcano_irfinder'

include { MULTIQC } from './modules/multiqc'

workflow {

    // 1. Build initial meta map channel from local FASTQ samplesheet
    ch_input_meta = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            def genome_info    = params.genomes ? params.genomes[params.genome] : null
            meta.gsm_id        = row.sample
            meta.sample_name   = row.sample
            meta.group1        = row.group
            meta.group2        = null
            meta.genome        = params.genome
            meta.star_index    = genome_info ? genome_info.star_index : null
            meta.bowtie2_index = genome_info ? genome_info.bowtie2_index : null
            
            def r1 = file(row.fastq_1)
            def r2 = row.fastq_2 ? file(row.fastq_2) : null
            meta.r1 = r1
            meta.r2 = r2
            meta.r3 = null
            meta.mode = r2 ? "PE" : "SE"

            return meta
        }

    // 2. Trimming
    ch_trimmed = TRIM_FASTQ(ch_input_meta)

    // 3. Update meta with trimmed file locations
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

    // 4. Alignments & Indexing (Separated)
    ch_star_te_out = STAR_TEALIGNMENT(ch_aligned_input)
    ch_indexed_bam = SAMTOOLS_INDEX(ch_star_te_out.bam)

    // 5. Downstream TE analyses
    ch_te_input = ch_indexed_bam.indexed_bam.map { meta, bam, bai ->
        def m = meta.clone()
        def selected_genome = params.genomes[meta.genome]
        m.bam         = bam
        m.gtf         = selected_genome.gtf
        m.te_gtf      = selected_genome.te_gtf
        m.te_loc      = selected_genome.te_loc
        m.scTE_idx    = selected_genome.scTE_idx
        m.scTE_tx_idx = selected_genome.scTE_tx_idx
        return m
    }

    // 6. Run TECOUNT analysis and summary
    ch_tecount_out = TECOUNT(ch_te_input)

    ch_all_counts = ch_tecount_out
        .map { meta, count_file -> count_file }
        .collect()

    ch_merged_matrix = MERGE_TECOUNTS(ch_all_counts)

    DESEQ2_TECOUNT(
        ch_merged_matrix.matrix,
        file(params.contrast_sheet)
    )

    ch_te_results_flat = DESEQ2_TECOUNT.out.results.flatten()

    ch_te_volcano_input = ch_te_results_flat
        .map { file ->
            def name = file.name.replaceAll("_deseq2_results\\.csv", "")
            return tuple(name, file)
        }

    VOLCANO_TECOUNT(ch_te_volcano_input)

    HEATMAP_TECOUNT(
        DESEQ2_TECOUNT.out.normalized_counts,
        file(params.contrast_sheet),
        DESEQ2_TECOUNT.out.combined_results
    )

    // Run coding gene analysis and summary based on TECOUNT
    DESEQ2_CODING(
        ch_merged_matrix.matrix,
        file(params.contrast_sheet)
    )

    ch_coding_results_flat = DESEQ2_CODING.out.results.flatten()

    ch_volcano_input = ch_coding_results_flat
        .map { file ->
            def name = file.name.replaceAll("_coding_deseq2_results\\.csv", "")
            return tuple(name, file)
        }

    VOLCANO_PLOT(ch_volcano_input)

    PATHWAY_ANALYSIS(ch_coding_results_flat)

    HEATMAP_ANALYSIS(
        DESEQ2_CODING.out.normalized_counts,
        file(params.contrast_sheet),
        DESEQ2_CODING.out.combined_results
    )

    // 7. Run TELOCAL process and downstream analysis
    ch_telocal_out = TELOCAL(ch_te_input)

    ch_all_telocal_counts = ch_telocal_out
        .map { meta, count_file -> count_file }
        .collect()

    ch_merged_telocal_matrix = MERGE_TELOCAL(ch_all_telocal_counts)

    EDGER_TELOCAL(ch_merged_telocal_matrix.matrix, file(params.contrast_sheet))

    // 8. Run scTE / scTELocal analyses
    SC_TE(ch_te_input)
    SC_TELOCAL(ch_te_input)

    // 9. Build inputs for IRFinder using trimmed FASTQ paths
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

    ch_irfinder_out = IRFINDER_FASTQ(ch_irfinder_input)
    ch_all_ir_dirs = ch_irfinder_out.ir_dir.collect()
    ch_merged_ir = MERGE_IRFINDER(ch_all_ir_dirs)
    DESEQ2_IRFINDER(ch_merged_ir.intron_matrix, ch_merged_ir.splice_matrix, file(params.contrast_sheet))
    
    HEATMAP_IRFINDER(
        DESEQ2_IRFINDER.out.ratio_matrix,
        file(params.contrast_sheet),
        DESEQ2_IRFINDER.out.combined_results
    )

    ch_ir_results_flat = DESEQ2_IRFINDER.out.results.flatten()
    ch_ir_volcano_input = ch_ir_results_flat
        .map { file ->
            def name = file.name.replaceAll("_irfinder_deseq2_results\\.csv", "")
            return tuple(name, file)
        }
    VOLCANO_IRFINDER(ch_ir_volcano_input)

    // 10. Collect QC logs for MultiQC
    ch_trim_qc = ch_trimmed.map { meta, files -> 
        (files instanceof List ? files : [files]).collect { it.getParent() }
    }.flatten()

    ch_star_qc = ch_star_te_out.bam.map { meta, bam -> 
        bam.getParent() 
    }

    ch_multiqc_inputs = ch_trim_qc
        .mix(ch_star_qc)
        .unique()
        .collect()

    MULTIQC(ch_multiqc_inputs)
}