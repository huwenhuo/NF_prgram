process TECOUNT {
    tag { meta.gsm_id }
    cpus 2
    memory '20 GB'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.tecount.cntTable")

    script:
    """
    TEcount \
        --sortByPos --format BAM --mode multi \
        -b ${meta.bam} \
        --GTF ${meta.gtf} \
        --TE ${meta.te_gtf} \
        --project ${meta.gsm_id}.tecount
    """
}
