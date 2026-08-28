process TELOCAL {
    tag { meta.gsm_id }
    cpus 1
    memory '10 GB'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.telocal.cntTable")

    script:
    """
    TElocal \
        --sortByPos -b ${meta.bam} \
        --GTF ${meta.gtf} \
        --TE ${meta.te_loc} \
        --stranded reverse \
        --project ${meta.gsm_id}.telocal
    """
}
