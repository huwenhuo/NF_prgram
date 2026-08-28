process DESEQ2_TECOUNT {
    tag "contrast_analysis_te"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/tecount_analysis", mode: 'copy'

    input:
    path counts_matrix
    path contrast_sheet

    output:
    path "*_deseq2_results.csv"         , emit: results
    path "all_deseq2_results.csv"       , emit: combined_results
    path "*_pca.pdf"                    , emit: pca_plots
    path "tecount_normalized_counts.csv", emit: normalized_counts

    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(ggplot2)
    library(data.table)

    meta_df <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    counts  <- read.table("${counts_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    is_ensembl_gene <- grepl("^ENSG|^ENSMUSG", rownames(counts), ignore.case = TRUE)
    counts <- counts[!is_ensembl_gene, , drop = FALSE]

    # Export full normalized matrix across valid samples
    valid_all_samples <- intersect(meta_df\$gsm_id, colnames(counts))
    full_counts <- counts[, valid_all_samples, drop = FALSE]
    full_meta <- meta_df[meta_df\$gsm_id %in% valid_all_samples, , drop = FALSE]
    rownames(full_meta) <- full_meta\$gsm_id

    dds_full <- DESeqDataSetFromMatrix(countData = round(full_counts), colData = full_meta, design = ~ 1)
    dds_full <- estimateSizeFactors(dds_full)
    write.csv(counts(dds_full, normalized = TRUE), file = "tecount_normalized_counts.csv", quote = FALSE)

    ctrl_patterns <- c("wt", "ctrl", "control")
    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file")
    candidate_cols <- setdiff(colnames(meta_df), reserved_cols)
    candidate_cols <- candidate_cols[grepl("^group|condition|treatment", candidate_cols, ignore.case = TRUE)]

    all_results_list <- list()

    for (col_name in candidate_cols) {
        sub_meta <- meta_df[!is.na(meta_df[[col_name]]) & tolower(as.character(meta_df[[col_name]])) != "na" & meta_df[[col_name]] != "", ]
        valid_samples <- intersect(sub_meta\$gsm_id, colnames(counts))
        if (length(valid_samples) < 2) next

        rownames(sub_meta) <- sub_meta\$gsm_id
        sub_meta <- sub_meta[valid_samples, , drop = FALSE]
        sub_counts <- counts[, valid_samples, drop = FALSE]

        raw_values <- as.character(sub_meta[[col_name]])
        sub_meta\$target_factor <- as.factor(ifelse(tolower(raw_values) %in% ctrl_patterns, "ctrl", raw_values))
        if (!"ctrl" %in% levels(sub_meta\$target_factor)) next

        sub_meta\$target_factor <- relevel(sub_meta\$target_factor, ref = "ctrl")
        treatments <- setdiff(levels(sub_meta\$target_factor), "ctrl")
        if (length(treatments) == 0) next

        dds <- DESeqDataSetFromMatrix(countData = round(sub_counts), colData = sub_meta, design = ~ target_factor)
        dds <- DESeq(dds)

        pdf(paste0(col_name, "_pca.pdf"))
        vsd <- vst(dds, blind = FALSE)
        print(plotPCA(vsd, intgroup = "target_factor") + ggtitle(paste("PCA (TE):", col_name)))
        dev.off()

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

    if (length(all_results_list) > 0) {
        fwrite(rbindlist(all_results_list, fill = TRUE), file = "all_deseq2_results.csv")
    }
    """
}