process GENERATE_BIGWIG {
    tag { meta.gsm_id }
    cpus 4
    memory 16.GB

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.bw")

    script:
    """
    bamCoverage \
        --bam ${meta.filt_bam} \
        --outFileName ${meta.gsm_id}.bw \
        --outFileFormat bigwig \
        --numberOfProcessors ${task.cpus} \
        --normalizeUsing CPM \
        --binSize 10 
    """
}
