nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet.tsv"  
params.fastq_dir = "/archive/InternalMedicine/Chung_lab/shared/sc/dna_meth/25075-03-10022025_141321"
params.outdir = "results"

selected_genome = params.genomes[params.genome]

//println "STAR index: ${selected_genome.star_index}"
//println "FASTA: ${selected_genome.fasta}"
//println "GTF: ${selected_genome.gtf}"

workflow { 

    // trimmed_fastq_ch should be produced by your TRIM_FASTQ
    // Example: trimmed_fastq_ch = samples_ch | TRIM_FASTQ
    trimmed_fastq_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def sample_id = row.Sample_Name
            def fq1 = file("${params.fastq_dir}/${sample_id}_*_R1_001.fastq.gz").first()
            def fq2 = file("${params.fastq_dir}/${sample_id}_*_R2_001.fastq.gz").first()
            tuple(sample_id, fq1, fq2)
        }
        | TRIM_FASTQ

    // Alignment and downstream
    cutrun_bam_ch = trimmed_fastq_ch | CUTRUN_ALIGN | SORT_INDEX | MARKDUP 

    // Peak calling and bigWig
    cutrun_bam_ch
        | CALL_PEAKS
        | PUBLISH_RESULTS

}

process TRIM_FASTQ {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.trimmed.fastq.gz"),
          path("${sample_id}_R2.trimmed.fastq.gz")

    script:
    """
    fastp \\
        -i ${R1} \\
        -I ${R2} \\
        -o ${sample_id}_R1.trimmed.fastq.gz \\
        -O ${sample_id}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --length_required 36 \\
        --thread 4 \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json
    """
}

process CUTRUN_ALIGN {
    tag { sample_id }
    cpus 8
    memory 32.GB

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam")

    // publish per-sample alignments optionally (controlled later)
    script:
    """
    set -euo pipefail

    # load module if needed (uncomment in cluster)
    # module load bowtie2/2.3.5

    # ensure bowtie2 index is set in selected_genome:
    if [ -z "${selected_genome.bowtie2_index}" ]; then
        echo "ERROR: selected_genome.bowtie2_index is not set" 1>&2
        exit 2
    fi

    # run bowtie2 (sensitive-local or --very-sensitive-local recommended for CUT&RUN)
    module load bowtie2/2.4.2
    module load samtools/1.6
    module load hisat2/2.1.0-intel
    bowtie2 \
        --threads ${task.cpus} \
        --very-sensitive-local \
        -x ${selected_genome.bowtie2_index} \
        -1 ${R1} -2 ${R2} \
        2> ${sample_id}_bowtie2.log | \
      samtools view -b -@ ${task.cpus} -o ${sample_id}.aligned.bam -
    """
}

                                                                                                                                                                                                        53,1          28%
process SORT_INDEX {
    tag { sample_id }
    cpus 4
    memory 16.GB

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam")

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.aligned.sorted.bam ${bam}
    samtools index ${sample_id}.aligned.sorted.bam
    """
}

process MARKDUP {
    tag { sample_id }
    cpus 4
    memory 16.GB

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.markdup.bam"), path("${sample_id}.aligned.sorted.markdup.bam.bai")

    script:
    """
    # Use samtools markdup (requires name-sorted then markdup) - perform name sort then markdup then coordinate sort
    module load samtools/1.6
    module load hisat2/2.1.0-intel
    samtools view -@ ${task.cpus} -F 0x4 -b ${sorted_bam} -o ${sample_id}.coord.bam

    # create name-sorted intermediate
    samtools sort -n -@ ${task.cpus} -o ${sample_id}.namesorted.bam ${sample_id}.coord.bam

    # mark duplicates
    samtools fixmate -m ${sample_id}.namesorted.bam ${sample_id}.fixmate.bam
    samtools sort -@ ${task.cpus} -o ${sample_id}.fixmate.sorted.bam ${sample_id}.fixmate.bam
    samtools markdup -r ${sample_id}.fixmate.sorted.bam ${sample_id}.aligned.sorted.markdup.bam

    # index
    samtools index ${sample_id}.aligned.sorted.markdup.bam
    """
}

process FILTER_BAM {
    tag { sample_id }
    cpus 2
    memory 8.GB

    input:
    tuple val(sample_id), path(markdup_bam)

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai")

    script:
    """

    module load samtools/1.6
    module load hisat2/2.1.0-intel

    # filter: remove duplicates (already removed with -r), keep MAPQ>=10 (or 30), remove chrM
    samtools view -@ ${task.cpus} -b -q 10 ${markdup_bam} \
        | samtools idxstats -  > /dev/null 2>&1

    # remove mitochondrial reads (common for CUT&RUN)
    samtools view -@ ${task.cpus} -b ${markdup_bam} \
        | samtools view -@ ${task.cpus} -h - \
        | awk 'BEGIN{OFS="\\t"} /^@/ {print; next} \$3!="chrM" && \$3!="MT" {print}' \
        | samtools view -@ ${task.cpus} -b -o ${sample_id}.ori.final.bam -

    samtools sort -@ ${task.cpus} -o ${sample_id}.final.bam ${sample_id}.ori.final.bam
    samtools index ${sample_id}.final.bam
    rm ${sample_id}.ori.final.bam
    """
}

process CALL_PEAKS {
    tag { sample_id }
    cpus 8
    memory 32.GB

    input:
    tuple val(sample_id), path(final_bam), path(final_bai)

    output:
    tuple val(sample_id), path("${sample_id}_peaks.narrowPeak"), path("${sample_id}.bigwig")

    publishDir "${params.outdir}/", mode: 'copy'

    script:
    """
    module load deeptools/3.5.0
    module load macs/2.1.2

    # create fragment BED for peak caller if needed
    # call MACS2 (BAMPE for paired-end CUT&RUN)
    macs2 callpeak -t ${final_bam} -f BAMPE -g ${selected_genome.genomeSize ?: 'mm'} -n ${sample_id} --outdir ./macs2_out --keep-dup all

    # macs2 will create *_peaks.narrowPeak in macs2_out
    cp macs2_out/${sample_id}_peaks.narrowPeak ${sample_id}_peaks.narrowPeak

    # create bigWig coverage using bamCoverage (deeptools)
    # normalize to CPM
    bamCoverage -b ${final_bam} -o ${sample_id}.bigwig --normalizeUsing CPM --binSize 10 -p ${task.cpus}
    """
}

process PUBLISH_RESULTS {
    tag { sample_id }
    input:
    tuple val(sample_id), path(peaks), path(bigwig)

    // final publish is handled in CALL_PEAKS already for sample-level outputs;
    // This step can be used to aggregate or further process multiple outputs.
    output:
    path "${sample_id}_peaks.narrowPeak", emit: peaks
    path "${sample_id}.bigwig", emit: bigwig

    script:
    """
    # noop: pass-through to keep outputs available to workflow
    """
}

                                                                                                                                                                                                        235,0-1       Bot
