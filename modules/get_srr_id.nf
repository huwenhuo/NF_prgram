process GET_SRR_ID {
    tag "${meta.gsm_id}"

    cpus 1
    memory '2 GB'

    input:
    val meta

    output:
    tuple val(meta), env(SRR_ID)

    script:
    """
    SRR_ID=\$(esearch -db sra -query "${meta.gsm_id}" | efetch -format runinfo | grep "SRR" | cut -d',' -f1 | head -n 1)
    """
}
