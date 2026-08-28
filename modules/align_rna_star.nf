process ALIGN_RNA_STAR {
    tag { meta.gsm_id }
    cpus 16
    memory 64.GB
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.Aligned.sortedByCoord.out.bam"), path("${meta.gsm_id}.ReadsPerGene.out.tab")

    script:
    def read_input = meta.trim_r2 ? "${meta.trim_r1},${meta.trim_r2}" : "${meta.trim_r1}"
    
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${meta.star_index} \
         --readFilesIn ${read_input} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${meta.gsm_id}. \
         --quantMode GeneCounts \
         --outStd Log
    """
}
