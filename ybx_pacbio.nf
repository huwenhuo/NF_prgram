#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet2.tsv"
params.outdir = "results"

selected_genome = params.genomes[params.genome]

process DOWNLOAD_SRA {
    tag "${sample}"
    cpus 4

    input:
        tuple val(sample), val(srr)

    output:
        tuple val(sample), path("*.fastq")

    script:
    """
    fasterq-dump ${srr} --split-3 -O ./
    """
}

process ALIGN_FASTQ {
    tag "${sample}"
    cpus 12 
    memory '80 GB'

    input:
        tuple val(sample), path(fastq)

    output:
        tuple val(sample), path("${sample}.sorted.bam")

    script:
    """
    module load samtools/1.19.2
    minimap2 -ax splice -uf -k14 -t ${task.cpus} ${selected_genome.fasta} ${fastq} > ${sample}.sam
    samtools view -bS ${sample}.sam | samtools sort -@ ${task.cpus} -m 4G -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    """
}

process REPEAT_COVERAGE {
    tag "${sample}"
    cpus 2
    memory '8 GB'

    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(sample), path(sorted_bam)

    output:
        tuple val(sample), path("${sample}_rmsk_coverage.txt"), path("${sample}_repeats_covered.bed")

    script:
    """
    module load bedtools/2.31.1

    # Compute coverage of repeats by reads
    bedtools coverage  -a ${selected_genome.te_bed} -b ${sorted_bam} > ${sample}_rmsk_coverage.txt
    bedtools intersect -a ${selected_genome.te_bed} -b ${sorted_bam} -wa -u > ${sample}_repeats_covered.bed

    """
}

// SQANTI3 isoform QC
process SQANTI3_QC {
    tag "${sample}"
    cpus 15
    memory '80 GB'

    publishDir "${params.outdir}/", mode: 'copy'

    input:
        tuple val(sample), path(sorted_bam)

    output:
        tuple val(sample), path("sqanti3_${sample}/")

    script:
    """
    ${params.img_sqanti3} sqanti3_qc.py ${sorted_bam} ${selected_genome.gtf} --output sqanti3_${sample} -p ${task.cpus}
    """
}

// ---- Workflow ----
workflow {

    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .set { samples_ch }

    bam_ch = samples_ch
        | map { row -> tuple(row.sample, row.srr) }
        | DOWNLOAD_SRA 
        | ALIGN_FASTQ

    bam_ch
        | SQANTI3_QC

    bam_ch
        | REPEAT_COVERAGE

}
