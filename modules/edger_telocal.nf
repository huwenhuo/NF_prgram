process EDGER_TELOCAL {
    tag "contrast_analysis_telocal"
    cpus 4
    memory '16 GB'
    publishDir "${params.output}/telocal_analysis", mode: 'copy'

    input:
    path counts_matrix
    path contrast_sheet

    output:
    path "*_telocal_edger_results.csv"  , emit: results
    path "all_telocal_edger_results.csv", emit: combined_results
    path "*_telocal_mds.pdf"            , emit: mds_plots

    script:
    """
    #!/usr/bin/env Rscript

    library(edgeR)
    library(ggplot2)
    library(data.table)

    # 1. Load inputs
    meta_df <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    counts  <- read.table("${counts_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

    # Convert character columns to numeric and replace any coercion NAs with 0
    counts[] <- lapply(counts, function(x) {
        num <- suppressWarnings(as.numeric(as.character(x)))
        num[is.na(num)] <- 0
        return(num)
    })
    counts <- as.matrix(counts)

    # 2. Exclude ENSEMBL coding genes (ENSG for Human, ENSMUSG for Mouse)
    is_ensembl_gene <- grepl("^ENSG|^ENSMUSG", rownames(counts), ignore.case = TRUE)
    counts <- counts[!is_ensembl_gene, , drop = FALSE]

    if (nrow(counts) == 0) {
        stop("No TElocal features remaining after removing ENSEMBL coding genes.")
    }

    ctrl_patterns <- c("wt", "ctrl", "control")
    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file", "telocal_file")

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

        # 3. Create DGEList and normalize
        dge <- DGEList(counts = sub_counts, group = sub_meta\$target_factor)
        
        # Filter low-count features to boost efficiency
        keep <- filterByExpr(dge)
        dge <- dge[keep, , keep.lib.sizes = FALSE]

        dge <- calcNormFactors(dge)

        # Output MDS Plot (edgeR equivalent of PCA)
        pdf(paste0(col_name, "_telocal_mds.pdf"))
        plotMDS(dge, labels = sub_meta\$target_factor, main = paste("MDS:", col_name))
        dev.off()

        # 4. Fit quasi-likelihood GLM model (fast & robust for large count matrices)
        design <- model.matrix(~ target_factor, data = sub_meta)
        dge <- estimateDisp(dge, design)
        fit <- glmQLFit(dge, design)

        # Extract contrast results for each treatment group
        for (i in seq_along(treatments)) {
            trt_group <- treatments[i]
            
            # Target column index in model matrix
            coef_idx <- paste0("target_factortrt_", trt_group)
            if (!coef_idx %in% colnames(design)) {
                coef_idx <- i + 1
            }

            qlf <- glmQLFTest(fit, coef = coef_idx)
            res_table <- topTags(qlf, n = Inf)\$table

            # Convert to data.table
            res_dt <- as.data.table(res_table, keep.rownames = "telocal_id")

            # Standardize column names to match DESeq2 layout (log2FoldChange, pvalue, padj)
            setnames(res_dt, 
                     old = c("logFC", "PValue", "FDR"), 
                     new = c("log2FoldChange", "pvalue", "padj"), 
                     skip_absent = TRUE)

            res_dt[, trt := trt_group]
            res_dt[, contrast_column := col_name]

            # Reorder columns with telocal_id, trt, and contrast_column first
            setcolorder(res_dt, c("telocal_id", "trt", "contrast_column"))

            output_name <- paste0(col_name, "_", trt_group, "_vs_ctrl_telocal_edger_results.csv")
            fwrite(res_dt, file = output_name)

            all_results_list[[length(all_results_list) + 1]] <- res_dt
        }
    }

    # 5. Combine all results into a single data.table and export
    if (length(all_results_list) > 0) {
        combined_dt <- rbindlist(all_results_list, fill = TRUE)
        fwrite(combined_dt, file = "all_telocal_edger_results.csv")
    }
    """
}
