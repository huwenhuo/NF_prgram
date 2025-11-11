#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet_ybx_bulk.tsv"
params.fastq_dir = "/archive/InternalMedicine/Chung_lab/shared/sc/raw_fastq/"
params.outdir = "results"

selected_genome = params.genomes[params.genome]

process DOWNLOAD_SRA {
    tag "${sample}"
    cpus 4

    input:
        tuple val(sample), val(srr)

    output:
        tuple val(sample), path("${srr}_1.fastq"), path("${srr}_2.fastq")

    script:
    """
    fasterq-dump ${srr} --split-3 -O ./
    """
}

process TRIM_FASTQ {
    tag "${sample_id}"

    // resource
    cpus = 2
    memory = 10.GB
    time = 10.h

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq.gz"),
          path("${sample_id}_R2.trimmed.fastq.gz")

    script:
    """
    fastp \
        -i ${R1} \
        -I ${R2} \
        -o ${sample_id}_R1.trimmed.fastq.gz \
        -O ${sample_id}_R2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --length_required 36 \
        --thread ${task.cpus} \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}

process STAR_teALIGNMENT {
    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    // resource
    cpus = 10
    memory = 60.GB
    time = 10.h

    input:
    tuple val(sample_id), path(R1_fastq), path(R2_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), path("${sample_id}_Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    module load star/2.7.11b
    STAR --genomeDir ${selected_genome.star_index}  \
         --runThreadN ${task.cpus} \
         --runMode alignReads \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 1000 \
         --outSAMmultNmax -1 \
         --outFilterMismatchNoverLmax 0.06  \
         --outMultimapperOrder Random \
         --winAnchorMultimapNmax 1000 \
         --alignTranscriptsPerReadNmax 1000 \
         --alignMatesGapMax 350 \
         --readFilesIn ${R1_fastq} ${R2_fastq} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_
    module load samtools/1.19.2
    samtools index ${sample_id}_Aligned.sortedByCoord.out.bam
    """
}

process SC_TELOCAL {
    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    // resource
    cpus = 5
    memory = 80.GB
    time = 20.h

    input:
    tuple val(sample_id), path(bam), path(idx_file) 

    output:
    path "${sample_id}_scTE.csv", emit: scTE_dir

    script:
    """
    module load samtools/1.19.2
    scTE \
        -i ${bam} \
        -p ${task.cpus} \
        -x ${idx_file} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${sample_id}_scTEtx
    """
}

process SC_TE {
    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    // resource
    cpus = 5
    memory = 80.GB
    time = 20.h

    input:
    tuple val(sample_id), path(bam), path(idx_file) 

    output:
    path "${sample_id}_scTE.csv", emit: scTE_dir

    script:
    """
    module load samtools/1.19.2
    scTE \
        -i ${bam} \
        -p ${task.cpus} \
        -x ${idx_file} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${sample_id}_scTE
    """
}

workflow {

    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .set { samples_ch }
    
    star_ch = samples_ch
        | map { row -> tuple(row.sample, row.srr) }
        | DOWNLOAD_SRA
        | TRIM_FASTQ
        | STAR_teALIGNMENT

    star_ch
        | map { sample_id, bam, bai -> tuple(sample_id, bam, file(selected_genome.scTE_idx)) } 
        | SC_TE
    
    star_ch
        | map { sample_id, bam, bai -> tuple(sample_id, bam, file(selected_genome.scTE_tx_idx)) } 
        | SC_TELOCAL
}
