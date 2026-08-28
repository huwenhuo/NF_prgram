process DESEQ2_CODING {
    tag "contrast_analysis_coding"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/coding_gene/", mode: 'copy'

    input:
    path counts_matrix
    path contrast_sheet

    output:
    path "*_coding_deseq2_results.csv"  , emit: results
    path "all_coding_deseq2_results.csv", emit: combined_results
    path "*_coding_pca.pdf"             , emit: pca_plots
    path "deseq2_normalized_counts.csv" , emit: normalized_counts 

    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(ggplot2)
    library(data.table)
    library(AnnotationDbi)

    # 1. Load inputs
    meta_df <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    counts  <- read.table("${counts_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

    # 2. Retain ONLY ENSEMBL coding genes
    is_ensembl_gene <- grepl("^ENSG|^ENSMUSG", rownames(counts), ignore.case = TRUE)
    counts <- counts[is_ensembl_gene, , drop = FALSE]

    if (nrow(counts) == 0) {
        stop("No ENSEMBL coding genes (ENSG/ENSMUSG) found in the count matrix.")
    }

    # Match matrix columns to metadata gsm_id
    valid_all_samples <- intersect(meta_df\$gsm_id, colnames(counts))
    full_counts       <- counts[, valid_all_samples, drop = FALSE]
    full_meta         <- meta_df[meta_df\$gsm_id %in% valid_all_samples, , drop = FALSE]
    rownames(full_meta) <- full_meta\$gsm_id

    # 3. Compute normalized counts across ALL samples once
    dds_full <- DESeqDataSetFromMatrix(
        countData = round(full_counts),
        colData   = full_meta,
        design    = ~ 1
    )
    dds_full <- estimateSizeFactors(dds_full)
    norm_counts_all <- counts(dds_full, normalized = TRUE)
    write.csv(norm_counts_all, file = "deseq2_normalized_counts.csv", quote = FALSE)

    # 4. Dynamic OrgDb annotation detection
    sample_gene <- rownames(counts)[1]
    if (grepl("^ENSG", sample_gene, ignore.case = TRUE)) {
        library(org.Hs.eg.db)
        org_db <- org.Hs.eg.db
        message("Detected Human ENSEMBL IDs (ENSG). Using org.Hs.eg.db.")
    } else if (grepl("^ENSMUSG", sample_gene, ignore.case = TRUE)) {
        library(org.Mm.eg.db)
        org_db <- org.Mm.eg.db
        message("Detected Mouse ENSEMBL IDs (ENSMUSG). Using org.Mm.eg.db.")
    } else {
        org_db <- NULL
        message("Unrecognized ENSEMBL format. Skipping gene symbol mapping.")
    }

    ctrl_patterns <- c("wt", "ctrl", "control")
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

        valid_samples <- intersect(sub_meta\$gsm_id, colnames(counts))
        if (length(valid_samples) < 2) {
            message(paste("Skipping column:", col_name, "- insufficient samples matching gsm_id."))
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
            countData = round(sub_counts),
            colData   = sub_meta,
            design    = ~ target_factor
        )
        dds <- DESeq(dds)

        # PCA export
        pdf(paste0(col_name, "_coding_pca.pdf"))
        vsd <- vst(dds, blind = FALSE)
        print(plotPCA(vsd, intgroup = "target_factor") + ggtitle(paste("PCA (Coding Genes):", col_name)))
        dev.off()

        # Extract contrast tables and map gene symbols
        for (trt_group in treatments) {
            res <- results(dds, contrast = c("target_factor", trt_group, "ctrl"))
            res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")

            clean_ids <- sub("\\\\..*", "", res_dt\$gene_id)

            if (!is.null(org_db)) {
                res_dt[, symbol := mapIds(
                    org_db,
                    keys      = clean_ids,
                    column    = "SYMBOL",
                    keytype   = "ENSEMBL",
                    multiVals = "first"
                )]
            } else {
                res_dt[, symbol := NA_character_]
            }

            res_dt[, trt := trt_group]
            res_dt[, contrast_column := col_name]

            setcolorder(res_dt, c("gene_id", "symbol", "trt", "contrast_column"))

            output_name <- paste0(col_name, "_", trt_group, "_vs_ctrl_coding_deseq2_results.csv")
            fwrite(res_dt, file = output_name)

            all_results_list[[length(all_results_list) + 1]] <- res_dt
        }
    }

    # Save combined table
    if (length(all_results_list) > 0) {
        combined_dt <- rbindlist(all_results_list, fill = TRUE)
        fwrite(combined_dt, file = "all_coding_deseq2_results.csv")
    }
    """
}