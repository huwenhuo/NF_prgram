process MERGE_SCTELOCAL {
    tag "merge_scte_matrix"
    cpus 4
    memory '32 GB'
    publishDir "${params.output}/merged_matrix", mode: 'copy'

    input:
    path count_csvs, stageAs: "?/*"

    output:
    path "scte_merged_counts.txt", emit: matrix

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Find all staged scTE CSV outputs recursively
    files <- list.files(".", pattern = "_scTEtx\\\\.csv\$", recursive = TRUE, full.names = TRUE)

    if (length(files) == 0) {
        stop("No _scTEtx.csv files found for merging.")
    }

    merged_dt <- NULL

    for (f in files) {
        # Extract GSM ID from filename (e.g., GSM4912339_scTEtx.csv -> GSM4912339)
        gsm_id <- gsub("_scTEtx\\\\.csv\$", "", basename(f))
        
        # Read 2-column count table (feature ID and count)
        dt <- fread(f, header = TRUE)
        setnames(dt, 1:2, c("feature", gsm_id))

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
    fwrite(merged_dt, file = "scte_merged_counts.txt", sep = "\t", quote = FALSE)
    """
}
