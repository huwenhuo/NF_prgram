process MERGE_TELOCAL {
    tag "merge_telocal_matrix"
    cpus 2
    memory '16 GB'
    publishDir "${params.output}/telocal_analysis", mode: 'copy'

    input:
    path count_tables, stageAs: "?/*"

    output:
    path "telocal_merged_counts.txt", emit: matrix

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Find all staged telocal count files recursively
    files <- list.files(".", pattern = "\\\\.telocal\\\\.cntTable\$", recursive = TRUE, full.names = TRUE)

    if (length(files) == 0) {
        stop("No .telocal.cntTable files found for merging.")
    }

    merged_dt <- NULL

    for (f in files) {
        # Extract GSM ID from filename (e.g., GSM4912339.telocal.cntTable -> GSM4912339)
        gsm_id <- gsub("\\\\.telocal\\\\.cntTable\$", "", basename(f))
        
        # Read 2-column count table (feature ID and count)
        dt <- fread(f, header = FALSE)
        setnames(dt, c("feature", gsm_id))

        if (is.null(merged_dt)) {
            merged_dt <- dt
        } else {
            merged_dt <- merge(merged_dt, dt, by = "feature", all = TRUE)
        }
    }

    # Replace NA values with 0
    for (col in names(merged_dt)) {
        set(merged_dt, i = which(is.na(merged_dt[[col]])), j = col, value = 0)
    }

    # Save tab-delimited count matrix
    fwrite(merged_dt, file = "telocal_merged_counts.txt", sep = "\t", quote = FALSE)
    """
}
