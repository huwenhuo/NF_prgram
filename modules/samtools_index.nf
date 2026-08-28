process SAMTOOLS_INDEX {
    tag "${meta.gsm_id}"
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """
}
