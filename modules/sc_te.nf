process SC_TE {
    tag { meta.gsm_id }
    cpus 5
    memory '80 GB'

    input:
    val meta

    output:
    path "${meta.gsm_id}_scTE.csv", emit: scte_dir

    script:
    """
    /work/InternalMedicine/s184335/sc//repos/scTE/bin/scTE \
        -i ${meta.bam} \
        -p ${task.cpus} \
        -x ${meta.scTE_idx} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${meta.gsm_id}_scTE
    """
}
