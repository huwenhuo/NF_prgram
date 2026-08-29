process HEATMAP_TECOUNT {
    tag "heatmap_te_features"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/tecount_analysis", mode: 'copy'

    input:
    path norm_counts_file
    path contrast_sheet
    path deseq2_results

    output:
    path "*_te_features_heatmap.pdf"    , emit: pdf, optional: true
    path "sig_te_features_log2_matrix.csv", emit: matrix, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(ComplexHeatmap)
    library(circlize)

    norm_counts <- read.csv("${norm_counts_file}", row.names = 1, check.names = FALSE)

    detect_sep <- function(filepath) {
        first_line <- readLines(filepath, n = 1)
        if (grepl(",", first_line)) return(",")
        if (grepl("\t", first_line)) return("\t")
        return(",")
    }
    
    contrast_sep <- detect_sep("${contrast_sheet}")
    meta_df <- read.table("${contrast_sheet}", header = TRUE, sep = contrast_sep, stringsAsFactors = FALSE)

    results_dt  <- fread("${deseq2_results}")

    # Filter for significant TE features (padj < 0.05 & absolute log2FC > 1)
    sig_tes <- unique(results_dt[padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(padj) & !is.na(log2FoldChange), tecount])

    message(sprintf("Found %d unique TE features matching thresholds.", length(sig_tes)))

    if (length(sig_tes) == 0) {
        message("Skipping TE Heatmap: No features met filters.")
        q(save = "no", status = 0)
    }

    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file")
    candidate_cols <- setdiff(colnames(meta_df), reserved_cols)
    group_col <- candidate_cols[grepl("^group|condition|treatment", candidate_cols, ignore.case = TRUE)][1]

    sub_meta <- meta_df[!is.na(meta_df[[group_col]]) & meta_df[[group_col]] != "", ]
    valid_samples <- intersect(sub_meta\$gsm_id, colnames(norm_counts))
    
    rownames(sub_meta) <- sub_meta\$gsm_id
    sub_meta <- sub_meta[valid_samples, , drop = FALSE]
    sub_counts <- norm_counts[, valid_samples, drop = FALSE]

    valid_sig_tes <- intersect(sig_tes, rownames(sub_counts))
    sub_counts <- sub_counts[valid_sig_tes, , drop = FALSE]

    log2_mat <- log2(as.matrix(sub_counts) + 1)
    z_matrix <- t(scale(t(log2_mat)))
    z_matrix <- z_matrix[complete.cases(z_matrix), ]

    fwrite(as.data.table(z_matrix, keep.rownames = "tecount"), file = "sig_te_features_log2_matrix.csv")

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
        name = "Z-score (Log2 Norm)",
        col = col_fun,
        top_annotation = col_annotation,
        show_row_names = show_rows,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = paste0("TE Normalized Heatmap (FDR < 0.05, |log2FC| > 1, n=", nrow(z_matrix), ")"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8)
    )

    pdf("all_samples_sig_te_features_heatmap.pdf", width = 8, height = max(6, nrow(z_matrix) * 0.015))
    draw(ht)
    dev.off()
    """
}