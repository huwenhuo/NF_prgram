process MULTIQC {
    cpus 10
    memory 40.GB
    
    input:
    path qc_inputs

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}
