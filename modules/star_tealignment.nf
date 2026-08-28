process STAR_TEALIGNMENT {
    tag "${meta.gsm_id}"
    cpus 10
    memory '60 GB'
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}_TE_Aligned.sortedByCoord.out.bam"), emit: bam

    script:
    def read_input = meta.trim_r2 ? "${meta.trim_r1} ${meta.trim_r2}" : "${meta.trim_r1}"

    """
    STAR --genomeDir ${meta.star_index} \
         --runThreadN ${task.cpus} \
         --runMode alignReads \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 1000 \
         --outSAMmultNmax -1 \
         --outFilterMismatchNoverLmax 0.06 \
         --outMultimapperOrder Random \
         --winAnchorMultimapNmax 1000 \
         --alignTranscriptsPerReadNmax 1000 \
         --alignMatesGapMax 350 \
         --readFilesIn ${read_input} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${meta.gsm_id}_TE_
    """
}
