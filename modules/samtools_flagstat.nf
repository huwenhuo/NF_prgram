process SAMTOOLS_FLAGSTAT {
    tag { meta.gsm_id }
    cpus 1
    memory 2.GB

    input:
    tuple val(meta), path(bam)

    output:
    path "${meta.gsm_id}.flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${meta.gsm_id}.flagstat.txt
    """
}
