process BISMARK {
    tag { meta.gsm_id }
    cpus 4
    memory 70.GB
    module 'bismark/0.24.1'

    input:
    val meta 

    output:
    tuple val(meta),
      path("*.bismark.cov.gz"),
      path("*.bedGraph.gz"),
      path("*.html"),
      path("*.txt"),
      optional: true

    script:
    def read_input = meta.trim_r2 ? "-1 ${meta.trim_r1} -2 ${meta.trim_r2}" : "${meta.trim_r1}"
    def is_paired = meta.trim_r2 ? "--paired-end" : "--single-end"
    
    """
    # 1. Alignment (using multicore without the forbidden --basename)
    bismark \\
        --genome ${meta.bismark_index} \\
        --bowtie2 \\
        --multicore ${task.cpus} \\
        ${read_input}

    # Normalize Bismark's default output name to match meta.gsm_id
    mv *bismark_bt2*.bam ${meta.gsm_id}_bismark.bam

    # 2. Extract methylation using the normalized BAM name
    bismark_methylation_extractor \\
        ${is_paired} \\
        --bedGraph \\
        --gzip \\
        --multicore ${task.cpus} \\
        ${meta.gsm_id}_bismark.bam

    # 3. Generate HTML report for the sample
    bismark2report \
        --alignment_report *_report.txt \
        --splitting_report *_splitting_report.txt

    """
}