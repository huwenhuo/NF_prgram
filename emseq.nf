#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet.tsvbak"
params.fastq_dir = "/archive/InternalMedicine/Chung_lab/shared/sc/dna_meth/25075-03-10022025_141321"
params.outdir = "results"

selected_genome = params.genomes[params.genome]

//println "STAR index: ${selected_genome.star_index}"
//println "FASTA: ${selected_genome.fasta}"
//println "GTF: ${selected_genome.gtf}"


process TRIM_FASTQ {
    tag "${sample_id}"

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
        --thread 4 \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}

process ALIGN_BISULFITE {
    tag "${sample_id}"

    // resource
    cpus = 16
    memory = 100.GB
    time = 10.h

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("${sample_id}_bismark.bam")

    script:
    """
    module load bismark/0.21.0
    module load bowtie2/2.4.2 
    module load samtools/1.6
    module load hisat2/2.1.0-intel
    bismark \\
        --genome ${selected_genome.bismark_index}  \\
        -1 ${R1} -2 ${R2} \\
        -o ./ \\
        -N 1  \\
        -L 20 \\
        --multicore 4
    mv ${sample_id}_R1.trimmed_bismark_bt2_pe.bam ${sample_id}_bismark.bam
    samtools sort -n -o ${sample_id}_bismark.sorted.bam ${sample_id}_bismark.bam
    mv ${sample_id}_bismark.sorted.bam ${sample_id}_bismark.bam
    #samtools index ${sample_id}_bismark.bam
    """
}

process METHYLATION_EXTRACT {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}", mode: 'move'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path("*")

    script:
    """
    module load bismark/0.21.0
    bismark_methylation_extractor \
        --paired-end \
        --gzip \
        --bedGraph \
        --counts \
        --CX \
        ${bam_file}
    #mv *.CX_report.txt.gz ${sample_id}_CpG_report.txt.gz
    """
}

workflow {
     samples_ch = Channel
         .fromPath(params.samplesheet)
         .splitCsv(header:true, sep:'\t')
         .map { row ->
             def sample_id = row.Sample_Name
             def fq1 = file("${params.fastq_dir}/${sample_id}_*_R1_001.fastq.gz").first()
             def fq2 = file("${params.fastq_dir}/${sample_id}_*_R2_001.fastq.gz").first()

	     //println "Sample: $sample_id" 
	     //println "R1 path: ${fq1}" 
	     //println "R2 path: ${fq2}"

             tuple(sample_id, fq1, fq2)
         }
     
     trimmed_ch = samples_ch | TRIM_FASTQ
     bam_ch     = trimmed_ch | ALIGN_BISULFITE
     meth_ch    = bam_ch | METHYLATION_EXTRACT

}

