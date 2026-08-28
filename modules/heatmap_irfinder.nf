process HEATMAP_IRFINDER {
    tag "heatmap_irfinder"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/irfinder_analysis", mode: 'copy'

    input:
    path ratio_matrix_file
    path contrast_sheet
    path deseq2_results

    output:
    path "*_irfinder_features_heatmap.pdf", emit: pdf, optional: true
    path "sig_irfinder_ratio_matrix.csv"  , emit: matrix, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(ComplexHeatmap)
    library(circlize)

    ratio_mat  <- read.csv("${ratio_matrix_file}", row.names = 1, check.names = FALSE)
    meta_df    <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    results_dt <- fread("${deseq2_results}")

    # Filter significant IR events (padj < 0.05 & absolute log2FC > 1)
    sig_irs <- unique(results_dt[padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(padj), ir_id])

    message(sprintf("Found %d unique IRFinder features matching thresholds.", length(sig_irs)))

    if (length(sig_irs) == 0) {
        message("Skipping IRFinder Heatmap: No features met filters.")
        q(save = "no", status = 0)
    }

    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file", "telocal_file", "scte_file")
    candidate_cols <- setdiff(colnames(meta_df), reserved_cols)
    group_col <- candidate_cols[grepl("^group|condition|treatment", candidate_cols, ignore.case = TRUE)][1]

    sub_meta <- meta_df[!is.na(meta_df[[group_col]]) & meta_df[[group_col]] != "", ]
    valid_samples <- intersect(sub_meta\$gsm_id, colnames(ratio_mat))
    
    rownames(sub_meta) <- sub_meta\$gsm_id
    sub_meta <- sub_meta[valid_samples, , drop = FALSE]
    sub_ratio <- ratio_mat[, valid_samples, drop = FALSE]

    valid_sig_irs <- intersect(sig_irs, rownames(sub_ratio))
    sub_ratio <- sub_ratio[valid_sig_irs, , drop = FALSE]

    # Z-score scaling across samples for pattern visualization
    z_matrix <- t(scale(t(as.matrix(sub_ratio))))
    z_matrix <- z_matrix[complete.cases(z_matrix), ]

    fwrite(as.data.table(z_matrix, keep.rownames = "ir_id"), file = "sig_irfinder_ratio_matrix.csv")

    annotation_df <- data.frame(Group = sub_meta[[group_col]])
    rownames(annotation_df) <- rownames(sub_meta)
    
    col_annotation <- HeatmapAnnotation(
        df = annotation_df,
        col = list(Group = structure(rainbow(length(unique(annotation_df\$Group))), 
                                      names = unique(annotation_df\$Group)))
    )

    col_fun <- colorRamp2(c(-2, 0, 2), c("#377EB8", "#FFFFFF", "#E41A1C"))
    show_rows <- nrow(z_matrix) <= 50

    ht <- Heatmap(
        z_matrix,
        name = "Z-score (IR Ratio)",
        col = col_fun,
        top_annotation = col_annotation,
        show_row_names = show_rows,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = paste0("IRFinder Heatmap (FDR < 0.05, |log2FC| > 1, n=", nrow(z_matrix), ")"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8)
    )

    pdf("all_samples_sig_irfinder_features_heatmap.pdf", width = 8, height = max(6, nrow(z_matrix) * 0.015))
    draw(ht)
    dev.off()
    """
}