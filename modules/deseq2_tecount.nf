process DESEQ2_TECOUNT {
    tag "contrast_analysis"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/deseq2_tecount", mode: 'copy'

    input:
    path counts_matrix
    path contrast_sheet

    output:
    path "*_deseq2_results.csv"  , emit: results
    path "all_deseq2_results.csv", emit: combined_results
    path "*_pca.pdf"             , emit: pca_plots

    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(ggplot2)
    library(data.table)

    # 1. Load inputs
    meta_df <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    counts  <- read.table("${counts_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    is_ensembl_gene <- grepl("^ENSG|^ENSMUSG", rownames(counts), ignore.case = TRUE)
    counts <- counts[!is_ensembl_gene, , drop = FALSE]

    ctrl_patterns <- c("wt", "ctrl", "control")
    # 'gsm_id' is now treated as the primary identifier instead of a reserved exclusion column
    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file")

    all_cols <- colnames(meta_df)
    candidate_cols <- setdiff(all_cols, reserved_cols)
    candidate_cols <- candidate_cols[grepl("^group|condition|treatment", candidate_cols, ignore.case = TRUE)]

    all_results_list <- list()

    for (col_name in candidate_cols) {

        # Filter out NA and blank rows
        sub_meta <- meta_df[!is.na(meta_df[[col_name]]) & 
                            tolower(as.character(meta_df[[col_name]])) != "na" & 
                            meta_df[[col_name]] != "", ]

        # Match using gsm_id against matrix columns
        valid_samples <- intersect(sub_meta\$gsm_id, colnames(counts))
        if (length(valid_samples) < 2) {
            message(paste("Skipping column:", col_name, "- insufficient samples matching gsm_id."))
            message(paste("GSM IDs in sheet:", paste(sub_meta\$gsm_id, collapse = ", ")))
            message(paste("Matrix columns:", paste(colnames(counts), collapse = ", ")))
            next
        }

        rownames(sub_meta) <- sub_meta\$gsm_id
        sub_meta   <- sub_meta[valid_samples, , drop = FALSE]
        sub_counts <- counts[, valid_samples, drop = FALSE]

        # Standardize control aliases to 'ctrl'
        raw_values <- as.character(sub_meta[[col_name]])
        standardized_values <- ifelse(tolower(raw_values) %in% ctrl_patterns, "ctrl", raw_values)
        sub_meta\$target_factor <- as.factor(standardized_values)

        if (!"ctrl" %in% levels(sub_meta\$target_factor)) {
            message(paste("Skipping column:", col_name, "- no control group (wt/ctrl/control) found."))
            next
        }

        sub_meta\$target_factor <- relevel(sub_meta\$target_factor, ref = "ctrl")
        treatments <- setdiff(levels(sub_meta\$target_factor), "ctrl")
        
        if (length(treatments) == 0) {
            message(paste("Skipping column:", col_name, "- no treatment groups found."))
            next
        }

        # Fit model
        dds <- DESeqDataSetFromMatrix(
            countData = sub_counts,
            colData   = sub_meta,
            design    = ~ target_factor
        )
        dds <- DESeq(dds)

        # PCA export
        pdf(paste0(col_name, "_pca.pdf"))
        vsd <- vst(dds, blind = FALSE)
        print(plotPCA(vsd, intgroup = "target_factor") + ggtitle(paste("PCA:", col_name)))
        dev.off()

        # Extract contrast tables
        for (trt_group in treatments) {
            res <- results(dds, contrast = c("target_factor", trt_group, "ctrl"))
            
            res_dt <- as.data.table(as.data.frame(res), keep.rownames = "tecount")
            res_dt[, trt := trt_group]
            res_dt[, contrast_column := col_name]
            setcolorder(res_dt, c("tecount", "trt", "contrast_column"))

            output_name <- paste0(col_name, "_", trt_group, "_vs_ctrl_deseq2_results.csv")
            fwrite(res_dt, file = output_name)

            all_results_list[[length(all_results_list) + 1]] <- res_dt
        }
    }

    # Save combined table
    if (length(all_results_list) > 0) {
        combined_dt <- rbindlist(all_results_list, fill = TRUE)
        fwrite(combined_dt, file = "all_deseq2_results.csv")
    }
    """
}
