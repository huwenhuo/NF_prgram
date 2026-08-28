process ALIGN_DNA {
    tag { meta.gsm_id }
    cpus 8
    memory 32.GB

    input:
    val meta 

    output:
    tuple val(meta), path("${meta.gsm_id}.sorted.bam"), path("${meta.gsm_id}.sorted.bam.bai")

    script:
    // Determine if we are in PE or SE mode for Bowtie2
    def read_input = meta.trim_r2 ? "-1 ${meta.trim_r1} -2 ${meta.trim_r2}" : "-U ${meta.trim_r1}"
    
    """
    bowtie2 --threads ${task.cpus} \
            --very-sensitive-local \
            --rg-id ${meta.gsm_id} \
            --rg "SM:${meta.gsm_id}" \
            --rg "PL:ILLUMINA" \
            -x ${meta.bowtie2_index} \
            ${read_input} \
            2> ${meta.gsm_id}_bowtie2.log | \
    samtools sort -@ ${task.cpus} -o ${meta.gsm_id}.sorted.bam -

    samtools index ${meta.gsm_id}.sorted.bam
    """
}
