process IRFINDER_FASTQ {
    tag "${meta.gsm_id}"
    cpus 5
    memory '40 GB'

    input:
    val meta

    output:
    path "ir_out_${meta.gsm_id}",       emit: ir_dir
    path "ir_out_${meta.gsm_id}/*.txt", emit: ir_results

    script:
    def reads = meta.trim_r2 ? "${meta.trim_r1} ${meta.trim_r2}" : "${meta.trim_r1}"

    """
    IRFinder -m FASTQ \
        -r ${meta.irfinder_index} \
        -d ir_out_${meta.gsm_id} \
        -t ${task.cpus} \
        ${reads}
    """
}
