process SC_TELOCAL {
    tag { meta.gsm_id }
    cpus 5
    memory '120 GB'

    input:
    val meta

    output:
    path "${meta.gsm_id}_scTEtx.csv", emit: scte_dir

    script:
    """
    /work/InternalMedicine/s184335/sc//repos/scTE/bin/scTE \
        -i ${meta.bam} \
        -p ${task.cpus} \
        -x ${meta.scTE_tx_idx} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${meta.gsm_id}_scTEtx
    """
}
