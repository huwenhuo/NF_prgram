process MARK_DUPLICATES {
    tag { meta.gsm_id }
    cpus 2
    memory 10.GB

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.md.bam"), path("${meta.gsm_id}.md.bam.bai"), path("${meta.gsm_id}.metrics.txt")

    script:
    """
    # Use the environment variable provided by the module
    java -Xmx8G -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        I=${meta.bam} \
        O=${meta.gsm_id}.md.bam \
        M=${meta.gsm_id}.metrics.txt \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT

    # Ensure index naming is consistent
    if [ ! -f "${meta.gsm_id}.md.bam.bai" ]; then
        mv ${meta.gsm_id}.md.bai ${meta.gsm_id}.md.bam.bai
    fi
    """
}
