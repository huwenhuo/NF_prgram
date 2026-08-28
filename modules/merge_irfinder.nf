process MERGE_IRFINDER {
    tag "merge_irfinder"
    cpus 2
    memory '16 GB'
    publishDir "${params.output}/irfinder_analysis", mode: 'copy'

    input:
    path ir_dirs, stageAs: "?/*"

    output:
    path "irfinder_intron_depth.txt"  , emit: intron_matrix
    path "irfinder_splice_coverage.txt", emit: splice_matrix

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Find all IRFinder-IR-nondir.txt files staged across directories
    files <- list.files(".", pattern = "IRFinder-IR-nondir\\\\.txt\$", recursive = TRUE, full.names = TRUE)

    if (length(files) == 0) {
        stop("No IRFinder-IR-nondir.txt files found.")
    }

    intron_list <- list()
    splice_list <- list()

    for (f in files) {
        # Extract gsm_id from parent folder path (ir_out_GSM12345 -> GSM12345)
        dir_name <- basename(dirname(f))
        gsm_id   <- gsub("^ir_out_", "", dir_name)

        # Read IRFinder output without header so columns are consistently V1, V2, etc.
        dt <- fread(f, header = FALSE, fill = TRUE)

        # Drop header row if present
        if (dt[1, 1] == "Chr" || dt[1, 1] == "#Chr") {
            dt <- dt[-1]
        }

        # Build locus identifier: Gene/Chr:Start-End:Strand
        dt[, ir_id := paste0(V4, "/", V1, ":", V2, "-", V3, ":", V6)]

        # Extract Intron Depth (Col 9) and pmax of Splice Exon Left/Right (Cols 17 & 18)
        dt[, intron_depth := round(as.numeric(V9))]
        dt[, max_splice   := round(pmax(as.numeric(V17), as.numeric(V18), na.rm = TRUE))]

        # Store sample counts
        intron_dt <- dt[, .(ir_id, count = intron_depth)]
        setnames(intron_dt, "count", gsm_id)

        splice_dt <- dt[, .(ir_id, count = max_splice)]
        setnames(splice_dt, "count", gsm_id)

        intron_list[[gsm_id]] <- intron_dt
        splice_list[[gsm_id]] <- splice_dt
    }

    # Merge across all samples
    merged_intron <- Reduce(function(x, y) merge(x, y, by = "ir_id", all = TRUE), intron_list)
    merged_splice <- Reduce(function(x, y) merge(x, y, by = "ir_id", all = TRUE), splice_list)

    # Fill NA values with 0
    for (col in names(merged_intron)) {
        set(merged_intron, i = which(is.na(merged_intron[[col]])), j = col, value = 0)
        set(merged_splice, i = which(is.na(merged_splice[[col]])), j = col, value = 0)
    }

    # Export merged matrices
    fwrite(merged_intron, file = "irfinder_intron_depth.txt", sep = "\t", quote = FALSE)
    fwrite(merged_splice, file = "irfinder_splice_coverage.txt", sep = "\t", quote = FALSE)
    """
}
