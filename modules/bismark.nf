process BISMARK {
    tag { meta.gsm_id }
    cpus 8
    memory 32.GB
    module 'bismark/0.24.1'

    input:
    val meta 

    output:
    tuple val(meta), path("${meta.gsm_id}_bismark*.bam"), path("${meta.gsm_id}_bismark*.bai"), optional: true
    path "*.html"

    script:
    def read_input = meta.trim_r2 ? "-1 ${meta.trim_r1} -2 ${meta.trim_r2}" : "${meta.trim_r1}"
    
    """
    bismark \
        --genome ${meta.bismark_index} \
        --bowtie2 \
        --multicore ${task.cpus} \
        ${read_input}

    # Add sorting and indexing commands matching your output expectations here
    """
}