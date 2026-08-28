process HEATMAP_ANALYSIS {
    tag "heatmap_de_genes"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/coding_gene", mode: 'copy'

    input:
    path norm_counts_file
    path contrast_sheet
    path deseq2_results

    output:
    path "*_de_genes_heatmap.pdf"       , emit: pdf, optional: true
    path "sig_de_genes_log2_matrix.csv" , emit: matrix, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(AnnotationDbi)

    # 1. Read inputs
    norm_counts <- read.csv("${norm_counts_file}", row.names = 1, check.names = FALSE)
    meta_df     <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    results_dt  <- fread("${deseq2_results}")

    # Clean ENSEMBL IDs
    rownames(norm_counts) <- sub("\\\\..*", "", rownames(norm_counts))
    results_dt[, clean_gene_id := sub("\\\\..*", "", gene_id)]

    # 2. Extract significant DE genes (padj < 0.05 & absolute log2FC > 1)
    sig_genes <- unique(results_dt[padj < 0.05 & abs(log2FoldChange) > 1 & !is.na(padj) & !is.na(log2FoldChange), clean_gene_id])

    message(sprintf("Found %d unique genes with padj < 0.05 across contrasts.", length(sig_genes)))

    if (length(sig_genes) == 0) {
        message("Skipping Heatmap: No genes met the padj < 0.05 cutoff.")
        q(save = "no", status = 0)
    }

    # Match metadata samples to matrix columns
    reserved_cols  <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file")
    candidate_cols <- setdiff(colnames(meta_df), reserved_cols)
    group_col      <- candidate_cols[grepl("^group|condition|treatment", candidate_cols, ignore.case = TRUE)][1]

    sub_meta <- meta_df[!is.na(meta_df[[group_col]]) & meta_df[[group_col]] != "", ]
    valid_samples <- intersect(sub_meta\$gsm_id, colnames(norm_counts))

    rownames(sub_meta) <- sub_meta\$gsm_id
    sub_meta   <- sub_meta[valid_samples, , drop = FALSE]
    sub_counts <- norm_counts[, valid_samples, drop = FALSE]

    # Subset matrix for significant genes
    valid_sig_genes <- intersect(sig_genes, rownames(sub_counts))
    sub_counts      <- sub_counts[valid_sig_genes, , drop = FALSE]

    # 3. Log2 transform log2(norm_counts + 1) & Z-score scaling
    log2_mat <- log2(as.matrix(sub_counts) + 1)
    z_matrix <- t(scale(t(log2_mat)))
    z_matrix <- z_matrix[complete.cases(z_matrix), ]

    # 4. Map ENSEMBL IDs to Gene Symbols for Row Labels
    sample_gene <- rownames(z_matrix)[1]
    if (grepl("^ENSG", sample_gene, ignore.case = TRUE)) {
        library(org.Hs.eg.db)
        org_db <- org.Hs.eg.db
    } else if (grepl("^ENSMUSG", sample_gene, ignore.case = TRUE)) {
        library(org.Mm.eg.db)
        org_db <- org.Mm.eg.db
    } else {
        org_db <- NULL
    }

    if (!is.null(org_db)) {
        mapped_symbols <- mapIds(
            org_db,
            keys      = rownames(z_matrix),
            column    = "SYMBOL",
            keytype   = "ENSEMBL",
            multiVals = "first"
        )
        # Use Symbol if mapped, fallback to ENSEMBL ID if NA
        row_labels <- ifelse(!is.na(mapped_symbols), mapped_symbols, rownames(z_matrix))
        # Ensure row labels remain unique for heatmap rendering
        row_labels <- make.unique(row_labels)
    } else {
        row_labels <- rownames(z_matrix)
    }

    rownames(z_matrix) <- row_labels

    # Export normalized matrix
    fwrite(as.data.table(z_matrix, keep.rownames = "symbol"), file = "sig_de_genes_log2_matrix.csv")

    # 5. Build Heatmap Annotations & Render
    annotation_df <- data.frame(Group = sub_meta[[group_col]])
    rownames(annotation_df) <- rownames(sub_meta)

    col_annotation <- HeatmapAnnotation(
        df  = annotation_df,
        col = list(Group = structure(rainbow(length(unique(annotation_df\$Group))), 
                                      names = unique(annotation_df\$Group)))
    )

    col_fun   <- colorRamp2(c(-2, 0, 2), c("#377EB8", "#FFFFFF", "#E41A1C"))
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
        column_title = paste0("DESeq2 Normalized Heatmap (FDR < 0.05, n=", nrow(z_matrix), ")"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8)
    )

    pdf("all_samples_sig_de_genes_heatmap.pdf", width = 8, height = max(6, nrow(z_matrix) * 0.015))
    draw(ht)
    dev.off()
    """
}