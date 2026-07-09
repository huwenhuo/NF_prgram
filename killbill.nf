#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { DOWNLOAD_FASTQ } from './modules.nf'
include { TRIM_FASTQ } from './modules.nf'
include { ALIGN_RNA_STAR } from './modules.nf'
include { ALIGN_DNA } from './modules.nf'
include { MARK_DUPLICATES } from './modules.nf'
include { FILTER_BAM } from './modules.nf'
include { GENERATE_BIGWIG } from './modules.nf'
include { MACS3_CALLPEAK_NoCONTROL } from './modules.nf'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DNA; SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_RNA; MULTIQC } from './modules.nf'

workflow {

    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            [
                gsm_id: row.gsm_id,
                assay: row.assay,
                sample_name: row.sample_name,
                genome: params.genome,
                bowtie2_index: params.genomes[params.genome].bowtie2_index,
                star_index: params.genomes[params.genome].star_index
            ]
        }.set{ samples_ch }

    // ==========================================
    // # download fastq files
    // ==========================================
    DOWNLOAD_FASTQ(samples_ch)
        .map { meta, files -> 
            def new_meta = meta.clone() 
            
            new_meta.r1 = files.find { it.name.contains('_1') }
            new_meta.r2 = files.find { it.name.contains('_2') }
            new_meta.r3 = files.find { it.name.contains('_3') }
            new_meta.mode = (new_meta.r2) ? (new_meta.r3 ? "PE_plus_R3" : "PE") : "SE"
            
            return new_meta
        }
        .set { download_ch }

    // ==========================================
    // # trime PE fastq files
    // ==========================================
    TRIM_FASTQ(download_ch)
        .map { meta, reads ->
            def new_meta = meta.clone()
            
            if (reads instanceof List) {
                // Paired-end: reads is a list of two files
                new_meta.trim_r1 = reads[0]
                new_meta.trim_r2 = reads[1]
            } else {
                // Single-end: reads is just a single file
                new_meta.trim_r1 = reads
                new_meta.trim_r2 = null
            }
            
            return new_meta
        }
        .set { trim_ch }
    
    // ==========================================
    // # RNA Alignment
    // ==========================================
    trim_ch
        .filter { it.assay == 'RNA-Seq' } 
        .set { rna_samples_ch }
    
    ALIGN_RNA_STAR(rna_samples_ch)
        .map { meta, bam, counts -> 
            def new_meta = meta.clone()
            new_meta.bam = bam
            new_meta.counts = counts
            return new_meta
        }
        .set { aligned_rna_ch }

    // ==========================================
    // # DNA sample (Includes ChIP-Seq, ATAC-Seq, and 3C-seq)
    // ==========================================
    trim_ch
        .filter { it.assay == 'ChIP-Seq' || it.assay == '3C-Seq' || it.assay == 'ATAC-Seq' }
        .set { dna_align_ch }
        
    // ## Align DNA
    ALIGN_DNA(dna_align_ch)
        .map { meta, bam, bai -> 
            // Enrich meta with the resulting BAM paths for downstream use
            def new_meta = meta.clone()
            new_meta.bam = bam
            new_meta.bai = bai
            return new_meta
        }
        .set { align_dna_ch }

    // ## Mark Duplicates
    MARK_DUPLICATES(align_dna_ch)
        .map { meta, dedup_bam, dedup_bai, metrics -> 
            def new_meta = meta.clone()
            new_meta.dedup_bam = dedup_bam
            new_meta.dedup_bai = dedup_bai
            new_meta.dedup_metrics = metrics
            return new_meta
        }
        .set { dedup_ch }

    // ## Filter BAM
    FILTER_BAM(dedup_ch)
        .map { meta, filt_bam, filt_bai ->
            def new_meta = meta.clone()
            new_meta.filt_bam = filt_bam
            new_meta.filt_bai = filt_bai
            return new_meta
        }
        .set { filt_ch }

    // ## Generate BigWig
    GENERATE_BIGWIG(filt_ch)

    // ## Call Peaks (Excludes 3C-seq; runs on both ChIP-Seq and ATAC-Seq)
    filt_ch
        .filter { it.assay != '3C-Seq' }
        .set { peak_ch }

    MACS3_CALLPEAK_NoCONTROL(peak_ch)

    // ==========================================
    // # QUALITY CONTROL (QC) SECTION
    // ==========================================

    // Map to: [ meta, meta.filt_bam ]
    filt_ch
        .map { meta -> [ meta, meta.filt_bam ] }
        .set { dna_flagstat_in_ch }
    
    SAMTOOLS_FLAGSTAT_DNA(dna_flagstat_in_ch)

    // Map to: [ meta, meta.bam ]
    aligned_rna_ch
        .map { meta -> [ meta, meta.bam ] }
        .set { rna_flagstat_in_ch }

    SAMTOOLS_FLAGSTAT_RNA(rna_flagstat_in_ch)

    // Collect all files for MultiQC
    def qc_files_ch = dedup_ch.map { meta -> meta.dedup_metrics }
        .mix( SAMTOOLS_FLAGSTAT_DNA.out )
        .mix( SAMTOOLS_FLAGSTAT_RNA.out )
        .collect()

    // Generate the final report
    MULTIQC(qc_files_ch)
}