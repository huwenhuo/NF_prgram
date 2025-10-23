#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet.tsv"
params.fastq_dir = "/archive/InternalMedicine/Chung_lab/shared/sc/dna_meth/25075-03-10022025_141321"
params.outdir = "results"

selected_genome = params.genomes[params.genome]

//println "STAR index: ${selected_genome.star_index}"
//println "FASTA: ${selected_genome.fasta}"
//println "GTF: ${selected_genome.gtf}"


process TRIM_FASTQ {
    tag "${sample_id}"

    // resource
    cpus = 4
    memory = 20.GB
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
        --thread 4 \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}


process STAR_teALIGNMENT {
    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    // resource
    cpus = 16
    memory = 100.GB
    time = 10.h

    input:
    tuple val(sample_id), path(R1_fastq), path(R2_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    """
    module load star/2.7.10b
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
    """
}

process STAR_ALIGNMENT {
    tag "${sample_id}"
    //publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    // resource
    cpus = 16
    memory = 100.GB
    time = 10.h

    input:
    tuple val(sample_id), path(R1_fastq), path(R2_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_ReadsPerGene.out.tab")   

    script:
    """
    module load star/2.7.10b
    STAR --genomeDir ${selected_genome.star_index}  \
         --runThreadN ${task.cpus} \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 10 \
         --outFilterMismatchNoverLmax 0.06  \
         --quantMode GeneCounts \
         --readFilesIn ${R1_fastq} ${R2_fastq} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_
    """
}


process TECOUNT {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(GTF_FILE), path(TE_GTF_FILE)

    output:
    tuple val(sample_id), path("${sample_id}.tecount.cntTable")

    script:
    """
    module load singularity/default
    singularity exec ${params.img_tecount} TEcount \
        --sortByPos --format BAM --mode multi \
        -b ${bam_file} \
        --GTF ${GTF_FILE} \
        --TE ${TE_GTF_FILE} \
        --project ${sample_id}.tecount
    touch ${sample_id}.tecount.cntTable
    """
}


process TELOCAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(GTF_FILE), path(TE_loc_GTF_FILE)

    output:
    tuple val(sample_id), path("${sample_id}.telocal.cntTable")

    script:
    """
    module load singularity/default
    singularity exec ${params.img_telocal} TElocal \\
        --sortByPos -b ${bam_file} \\
        --GTF ${GTF_FILE} \\
        --TE ${TE_loc_GTF_FILE} \\
        --stranded reverse \\
        --project ${sample_id}.telocal
    touch ${sample_id}.telocal.cntTable 

    """
}

workflow {

    // ----------------------
    // Step 1: Trim FASTQs
    // ----------------------
    trimmed_fastq_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def sample_id = row.Sample_Name
            def fq1 = file("${params.fastq_dir}/${sample_id}_*_R1_001.fastq.gz")
            def fq2 = file("${params.fastq_dir}/${sample_id}_*_R2_001.fastq.gz")
            // Print for testing
            //println "Sample: $sample_id"
            //println "R1 path: ${fq1}"
            //println "R2 path: ${fq2}"
            

            tuple(sample_id, fq1, fq2)
        }
        | TRIM_FASTQ

    // ----------------------
    // Step 2a: STAR alignment for normal genes
    // ----------------------
    //star_gene_bam_ch = trimmed_fastq_ch
    //    | STAR_ALIGNMENT   // STAR run with --quantMode GeneCounts and low multimapping

    // ----------------------
    // Step 2b: STAR alignment for TE counting
    // ----------------------
    star_te_bam_ch = trimmed_fastq_ch
        | STAR_teALIGNMENT      // STAR run with --outFilterMultimapNmax 1000, no GeneCounts

    // ----------------------
    // Step 3: TEcount
    // ----------------------
    star_te_bam_ch
        .map { sample_id, bam ->
            tuple(sample_id, bam, file(selected_genome.gtf), file(selected_genome.te_gtf))
        }
        | TECOUNT

    // ----------------------
    // Step 3: TElocal
    // ----------------------
    //star_te_bam_ch
    //    .map { sample_id, bam ->
    //        tuple(sample_id, bam, file(selected_genome.gtf), file(selected_genome.te_loc))
    //    }
    //    | TELOCAL


}

