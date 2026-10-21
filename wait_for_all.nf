process WAIT_FOR_ALL {
    publishDir "${params.output}/pipeline_info", mode: 'copy', enabled: false
    
    input:
    val all_outputs

    output:
    val true, emit: complete

    script:
    """
    echo "All upstream sample processes completed successfully."
    """
}