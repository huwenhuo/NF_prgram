#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { DOWNLOAD_FASTQ } from './modules.nf'
include { TRIM_FASTQ_PE; TRIM_FASTQ_SE } from './modules.nf'
include { ALIGN_DNA_PE; ALIGN_DNA_SE } from './modules.nf'
include { ALIGN_RNA_STAR_PE; ALIGN_RNA_STAR_SE } from './modules.nf'
include { MARK_DUPLICATES } from './modules.nf'
include { FILTER_BAM } from './modules.nf'
include { GENERATE_BIGWIG } from './modules.nf'
include { MACS3_CALLPEAK_NoCONTROL } from './modules.nf'
include { TRIM_FASTQ_PE as TRIM_3C_PE; TRIM_FASTQ_SE as TRIM_3C_SE } from './modules.nf'
include { ALIGN_DNA_PE as ALIGN_3C_PE; ALIGN_DNA_SE as ALIGN_3C_SE } from './modules.nf'

workflow {
    // 1. Setup channel
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.gsm_id, row) }
        .set { samples_ch }

    // 2. Download
    DOWNLOAD_FASTQ(samples_ch.map { it[0] })

    // 3. Join and Branch by Assay AND File Count (PE vs SE)
    DOWNLOAD_FASTQ.out
        .join(samples_ch)
        .branch { gsm_id, files, row ->
            // files is a list, e.g., [file1.gz, file2.gz]
            println "DEBUG: Processing $gsm_id, File count: ${files.size()}, Files: $files"
            is_pe: files.size() > 1
            is_se: files.size() == 1
        }
        .set { data_ch }

    // 4. ChIP-Seq Path
    data_ch.is_pe.filter { it[2].assay == 'ChIP-Seq' }.set { chip_pe }
    data_ch.is_se.filter { it[2].assay == 'ChIP-Seq' }.set { chip_se }
    
    TRIM_FASTQ_PE(chip_pe.map { [it[0], it[1]] })
    TRIM_FASTQ_SE(chip_se.map { [it[0], it[1]] })
    
    ALIGN_DNA_PE(TRIM_FASTQ_PE.out, params.genomes[params.genome].bowtie2_index)
    ALIGN_DNA_SE(TRIM_FASTQ_SE.out, params.genomes[params.genome].bowtie2_index)
    
    // Combine back for downstream
    ALIGN_DNA_PE.out.mix(ALIGN_DNA_SE.out).set { bam_ch }
    
    MARK_DUPLICATES(bam_ch)
    FILTER_BAM(MARK_DUPLICATES.out)
    GENERATE_BIGWIG(FILTER_BAM.out)
    MACS3_CALLPEAK_NoCONTROL(FILTER_BAM.out)

    // 5. RNA-Seq Path
    data_ch.is_pe.filter { it[2].assay == 'RNA-Seq' }.set { rna_pe }
    data_ch.is_se.filter { it[2].assay == 'RNA-Seq' }.set { rna_se }
    
    ALIGN_RNA_STAR_PE(rna_pe.map { [it[0], it[1]] }, params.genomes[params.genome].star_index)
    ALIGN_RNA_STAR_SE(rna_se.map { [it[0], it[1]] }, params.genomes[params.genome].star_index)

    // 6. 3C-Seq Path
    data_ch.is_pe.filter { it[2].assay == '3C-Seq' }.set { c3_pe }
    data_ch.is_se.filter { it[2].assay == '3C-Seq' }.set { c3_se }
    
    TRIM_3C_PE(c3_pe.map { [it[0], it[1]] })
    TRIM_3C_SE(c3_se.map { [it[0], it[1]] })
    ALIGN_3C_PE(TRIM_3C_PE.out, params.genomes[params.genome].bowtie2_index)
    ALIGN_3C_SE(TRIM_3C_SE.out, params.genomes[params.genome].bowtie2_index)
}