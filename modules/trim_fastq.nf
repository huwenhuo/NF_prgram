process TRIM_FASTQ {
    tag { meta.gsm_id }
    cpus 4
    memory 8.GB

    input:
    val meta

    output:
    tuple val(meta), path("*.trimmed.fastq.gz")

    script:
    if (meta.mode == "PE") {
        """
        fastp \
            -i ${meta.r1} -I ${meta.r2} \
            -o ${meta.gsm_id}_R1.trimmed.fastq.gz -O ${meta.gsm_id}_R2.trimmed.fastq.gz \
            --disable_quality_filtering --length_required 20 --detect_adapter_for_pe \
            --thread ${task.cpus} \
            --html ${meta.gsm_id}_fastp.html --json ${meta.gsm_id}_fastp.json
        """
    } else {
        """
        fastp \
            -i ${meta.r1} \
            -o ${meta.gsm_id}_R1.trimmed.fastq.gz \
            --disable_quality_filtering --length_required 20 \
            --thread ${task.cpus} \
            --html ${meta.gsm_id}_fastp.html --json ${meta.gsm_id}_fastp.json
        """
    }
}
