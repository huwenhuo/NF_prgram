process DESEQ2_IRFINDER {
    tag "contrast_analysis_irfinder"
    cpus 4
    memory '32 GB'
    publishDir "${params.output}/deseq2_irfinder", mode: 'copy'

    input:
    path intron_matrix
    path splice_matrix
    path contrast_sheet

    output:
    path "*_irfinder_deseq2_results.csv"  , emit: results, optional: true
    path "all_irfinder_deseq2_results.csv", emit: combined_results, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(data.table)

    meta_df    <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    intron_mat <- read.table("${intron_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    splice_mat <- read.table("${splice_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

    # Ensure feature alignment
    common_ids <- intersect(rownames(intron_mat), rownames(splice_mat))
    intron_mat <- intron_mat[common_ids, , drop = FALSE]
    splice_mat <- splice_mat[common_ids, , drop = FALSE]

    ctrl_patterns <- c("wt", "ctrl", "control")
    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file", "telocal_file", "scte_file")

    all_cols <- colnames(meta_df)
    candidate_cols <- setdiff(all_cols, reserved_cols)
    candidate_cols <- candidate_cols[grepl("^group|condition|treatment", candidate_cols, ignore.case = TRUE)]

    all_results_list <- list()

    for (col_name in candidate_cols) {

        sub_meta <- meta_df[!is.na(meta_df[[col_name]]) & 
                            tolower(as.character(meta_df[[col_name]])) != "na" & 
                            meta_df[[col_name]] != "", ]

        valid_samples <- intersect(sub_meta\$gsm_id, colnames(intron_mat))

        if (length(valid_samples) < 2) next

        rownames(sub_meta) <- sub_meta\$gsm_id
        sub_meta <- sub_meta[valid_samples, , drop = FALSE]

        sub_intron <- intron_mat[, valid_samples, drop = FALSE]
        sub_splice <- splice_mat[, valid_samples, drop = FALSE]

        # Filter low counts
        keep <- rowSums(sub_intron) > 10
        sub_intron <- sub_intron[keep, , drop = FALSE]
        sub_splice <- sub_splice[keep, , drop = FALSE]

        if (nrow(sub_intron) == 0) next

        raw_values <- as.character(sub_meta[[col_name]])
        standardized_values <- ifelse(tolower(raw_values) %in% ctrl_patterns, "ctrl", raw_values)
        sub_meta\$target_factor <- as.factor(standardized_values)

        if (!"ctrl" %in% levels(sub_meta\$target_factor)) next

        sub_meta\$target_factor <- relevel(sub_meta\$target_factor, ref = "ctrl")
        treatments <- setdiff(levels(sub_meta\$target_factor), "ctrl")

        if (length(treatments) == 0) next

        # Prepare combined count matrix
        colnames(sub_intron) <- paste0("IR_", colnames(sub_intron))
        colnames(sub_splice) <- paste0("Splice_", colnames(sub_splice))
        combined_counts <- cbind(sub_intron, sub_splice)

        col_data <- data.frame(
            sampleID      = factor(rep(valid_samples, 2)),
            target_factor = factor(rep(as.character(sub_meta\$target_factor), 2)),
            IRFinder      = factor(c(rep("IR", length(valid_samples)), rep("Splice", length(valid_samples))), levels = c("Splice", "IR"))
        )
        rownames(col_data) <- colnames(combined_counts)

        # Standard non-paired IRFinder design formula
        dds <- DESeqDataSetFromMatrix(
            countData = combined_counts,
            colData   = col_data,
            design    = ~ target_factor + IRFinder + target_factor:IRFinder
        )

        sizeFactors(dds) <- rep(1, ncol(combined_counts))
        dds <- DESeq(dds, quiet = TRUE)

        # Extract contrast results for each treatment level vs ctrl
        for (trt_group in treatments) {
            target_coef <- paste0("target_factor", trt_group, ".IRFinderIR")
            
            if (!target_coef %in% resultsNames(dds)) {
                matched_coef <- grep(paste0("IRFinderIR.*", trt_group, "|", trt_group, ".*IRFinderIR"), resultsNames(dds), value = TRUE)
                if (length(matched_coef) > 0) {
                    target_coef <- matched_coef[1]
                } else {
                    next
                }
            }

            res <- results(dds, name = target_coef)
            res_dt <- as.data.table(as.data.frame(res), keep.rownames = "ir_id")
            
            # Robust split handling any number of slash-delimited components
            res_dt[, gene := sub("/.*", "", ir_id)]
            res_dt[, locus := sub("^[^/]*/", "", ir_id)]

            res_dt[, trt := trt_group]
            res_dt[, contrast_column := col_name]

            setcolorder(res_dt, c("ir_id", "gene", "locus", "trt", "contrast_column"))

            output_name <- paste0(col_name, "_", trt_group, "_vs_ctrl_irfinder_deseq2_results.csv")
            fwrite(res_dt, file = output_name)

            all_results_list[[length(all_results_list) + 1]] <- res_dt
        }
    }

    if (length(all_results_list) > 0) {
        combined_dt <- rbindlist(all_results_list, fill = TRUE)
        fwrite(combined_dt, file = "all_irfinder_deseq2_results.csv")
    }
    """
}
