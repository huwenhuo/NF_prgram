process FILTER_BAM {
    tag { meta.gsm_id }
    cpus 2
    memory 4.GB

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.filtered.bam"), path("${meta.gsm_id}.filtered.bam.bai")

    script:
    """
    # -q 30: Quality score >= 30
    # -F 4: Exclude unmapped reads
    # -F 256: Exclude secondary alignments
    # -F 1024: Exclude PCR/optical duplicates
    # -b: output BAM format
    
    samtools view -b -q 30 -F 1804 ${meta.dedup_bam} > ${meta.gsm_id}.filtered.bam
    
    # Index the filtered BAM
    samtools index ${meta.gsm_id}.filtered.bam
    """
}
