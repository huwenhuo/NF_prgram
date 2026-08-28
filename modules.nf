process GET_SRR_ID {
    tag "${meta.gsm_id}"

    cpus 1
    memory '2 GB'

    input:
    val meta

    output:
    tuple val(meta), env(SRR_ID)

    script:
    """
    SRR_ID=\$(esearch -db sra -query "${meta.gsm_id}" | efetch -format runinfo | grep "SRR" | cut -d',' -f1 | head -n 1)
    """
}

process DOWNLOAD_FASTQ {
    tag "${meta.gsm_id}"

    cpus 4
    memory '8 GB'

    input:
    tuple val(meta), val(srr_id)

    output:
    tuple val(meta), path("${meta.gsm_id}*_?.fastq.gz", arity: '1..*')

    script:
    """
    # 1. Prefetch the data 
    prefetch ${srr_id}
    
    # 2. Extract using fasterq-dump
    fasterq-dump --split-3 --include-technical --threads ${task.cpus} ${srr_id}

    # 3. Count how many fastq files were generated before renaming
    file_count=\$(ls ${srr_id}*.fastq 2>/dev/null | wc -l)

    # 4. Rename files from SRR ID to GSM ID
    for file in ${srr_id}*.fastq; do
        suffix=\${file#${srr_id}}
        mv "\$file" "${meta.gsm_id}\${suffix}"
    done
    
    # 5. Safety handling based on file count
    if [ "\$file_count" -eq 1 ]; then
        # True Single-End: rename gsm_id.fastq to gsm_id_1.fastq
        if [ -f "${meta.gsm_id}.fastq" ]; then
            mv "${meta.gsm_id}.fastq" "${meta.gsm_id}_1.fastq"
        fi
    elif [ -f "${meta.gsm_id}.fastq" ]; then
        # Paired-End with unmatched reads: isolate singletons
        mv "${meta.gsm_id}.fastq" "${meta.gsm_id}_unmatched.fastq"
    fi
    
    # 6. Gzip target biological reads
    gzip ${meta.gsm_id}_*.fastq
    """
}

process DOWNLOAD_FASTQ_old {
    tag "${meta.gsm_id}"

    cpus 2
    memory 8.GB
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}*_?.fastq.gz", arity: '1..*')

    script:
    """
    # 1. Retrieve the SRR ID
    SRR_ID=\$(esearch -db sra -query "${meta.gsm_id}" | efetch -format runinfo | grep "SRR" | cut -d',' -f1 | head -n 1)
    
    # 2. Prefetch the data 
    prefetch \$SRR_ID
    
    # 3. Extract using fasterq-dump
    fasterq-dump --split-3 --include-technical --threads 4 \$SRR_ID

    # 4. Count how many fastq files were generated before renaming
    file_count=\$(ls \${SRR_ID}*.fastq 2>/dev/null | wc -l)

    # 5. Rename files from SRR ID to GSM ID
    for file in \${SRR_ID}*.fastq; do
        suffix=\${file#\$SRR_ID}
        mv "\$file" "${meta.gsm_id}\${suffix}"
    done
    
    # 6. Safety handling based on file count
    if [ "\$file_count" -eq 1 ]; then
        # True Single-End: rename gsm_id.fastq to gsm_id_1.fastq
        if [ -f "${meta.gsm_id}.fastq" ]; then
            mv "${meta.gsm_id}.fastq" "${meta.gsm_id}_1.fastq"
        fi
    elif [ -f "${meta.gsm_id}.fastq" ]; then
        # Paired-End with unmatched reads: isolate the singletons so they don't break downstream patterns
        mv "${meta.gsm_id}.fastq" "${meta.gsm_id}_unmatched.fastq"
    fi
    
    # 7. Gzip only the biological reads we want to capture (_1, _2, _3)
    gzip ${meta.gsm_id}_*.fastq
    """
}

process TRIM_FASTQ {
    tag { meta.gsm_id }
    cpus 4
    memory 8.GB

    input:
    val meta

    output:
    tuple val(meta), path("*.trimmed.fastq.gz")

    script:
    if (meta.mode == "PE") {
        """
        fastp \
            -i ${meta.r1} -I ${meta.r2} \
            -o ${meta.gsm_id}_R1.trimmed.fastq.gz -O ${meta.gsm_id}_R2.trimmed.fastq.gz \
            --disable_quality_filtering --length_required 20 --detect_adapter_for_pe \
            --thread ${task.cpus} \
            --html ${meta.gsm_id}_fastp.html --json ${meta.gsm_id}_fastp.json
        """
    } else {
        """
        fastp \
            -i ${meta.r1} \
            -o ${meta.gsm_id}_R1.trimmed.fastq.gz \
            --disable_quality_filtering --length_required 20 \
            --thread ${task.cpus} \
            --html ${meta.gsm_id}_fastp.html --json ${meta.gsm_id}_fastp.json
        """
    }
}

process ALIGN_RNA_STAR {
    tag { meta.gsm_id }
    cpus 16
    memory 64.GB
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.Aligned.sortedByCoord.out.bam"), path("${meta.gsm_id}.ReadsPerGene.out.tab")

    script:
    def read_input = meta.trim_r2 ? "${meta.trim_r1},${meta.trim_r2}" : "${meta.trim_r1}"
    
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${meta.star_index} \
         --readFilesIn ${read_input} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${meta.gsm_id}. \
         --quantMode GeneCounts \
         --outStd Log
    """
}

process STAR_TEALIGNMENT {
    tag "${meta.gsm_id}"
    cpus 10
    memory '60 GB'
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}_TE_Aligned.sortedByCoord.out.bam"), emit: bam

    script:
    def read_input = meta.trim_r2 ? "${meta.trim_r1} ${meta.trim_r2}" : "${meta.trim_r1}"

    """
    STAR --genomeDir ${meta.star_index} \
         --runThreadN ${task.cpus} \
         --runMode alignReads \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 1000 \
         --outSAMmultNmax -1 \
         --outFilterMismatchNoverLmax 0.06 \
         --outMultimapperOrder Random \
         --winAnchorMultimapNmax 1000 \
         --alignTranscriptsPerReadNmax 1000 \
         --alignMatesGapMax 350 \
         --readFilesIn ${read_input} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${meta.gsm_id}_TE_
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.gsm_id}"
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index ${bam}
    """
}

process TECOUNT {
    tag { meta.gsm_id }
    cpus 2
    memory '20 GB'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.tecount.cntTable")

    script:
    """
    TEcount \
        --sortByPos --format BAM --mode multi \
        -b ${meta.bam} \
        --GTF ${meta.gtf} \
        --TE ${meta.te_gtf} \
        --project ${meta.gsm_id}.tecount
    """
}

process MERGE_TECOUNTS {
    tag "all_samples_tecount"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/tecount_matrix", mode: 'copy'

    input:
    path count_files

    output:
    path "raw_tecounts_matrix.tsv", emit: matrix

    script:
    """
    #!/usr/bin/env Rscript

    # 1. Use character class [.] to avoid backslash escaping issues in Nextflow
    files <- list.files(pattern = "[.]tecount", full.names = TRUE)

    # 2. Read each file and extract sample name from file prefix
    counts_list <- lapply(files, function(f) {
        sample_name <- sub("[.]tecount.*", "", basename(f))
        df <- read.table(f, header = TRUE, row.names = 1)
        colnames(df) <- sample_name
        return(df)
    })

    # 3. Merge all count columns into a single matrix
    merged_matrix <- do.call(cbind, counts_list)

    # 4. Save consolidated counts matrix
    write.table(merged_matrix, file = "raw_tecounts_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
    """
}

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

process DESEQ2_CODING {
    tag "contrast_analysis_coding"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/deseq2_coding", mode: 'copy'

    input:
    path counts_matrix
    path contrast_sheet

    output:
    path "*_coding_deseq2_results.csv"  , emit: results
    path "all_coding_deseq2_results.csv", emit: combined_results
    path "*_coding_pca.pdf"             , emit: pca_plots

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

    # 2. Retain ONLY ENSEMBL coding genes (ENSG for Human, ENSMUSG for Mouse)
    is_ensembl_gene <- grepl("^ENSG|^ENSMUSG", rownames(counts), ignore.case = TRUE)
    counts <- counts[is_ensembl_gene, , drop = FALSE]

    if (nrow(counts) == 0) {
        stop("No ENSEMBL coding genes (ENSG/ENSMUSG) found in the count matrix.")
    }

    # 3. Dynamic OrgDb annotation detection
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

        # Fit model
        dds <- DESeqDataSetFromMatrix(
            countData = sub_counts,
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

            # Clean ENSEMBL IDs (escaped backslashes for Nextflow template evaluation)
            clean_ids <- sub("\\\\..*", "", res_dt\$gene_id)

            # Map symbols using OrgDb package
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

            # Reorder columns with gene_id and symbol at the beginning
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

process TELOCAL {
    tag { meta.gsm_id }
    cpus 1
    memory '10 GB'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.telocal.cntTable")

    script:
    """
    TElocal \
        --sortByPos -b ${meta.bam} \
        --GTF ${meta.gtf} \
        --TE ${meta.te_loc} \
        --stranded reverse \
        --project ${meta.gsm_id}.telocal
    """
}

process MERGE_TELOCAL {
    tag "merge_telocal_matrix"
    cpus 2
    memory '16 GB'
    publishDir "${params.output}/merged_matrix", mode: 'copy'

    input:
    path count_tables, stageAs: "?/*"

    output:
    path "telocal_merged_counts.txt", emit: matrix

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Find all staged telocal count files recursively
    files <- list.files(".", pattern = "\\\\.telocal\\\\.cntTable\$", recursive = TRUE, full.names = TRUE)

    if (length(files) == 0) {
        stop("No .telocal.cntTable files found for merging.")
    }

    merged_dt <- NULL

    for (f in files) {
        # Extract GSM ID from filename (e.g., GSM4912339.telocal.cntTable -> GSM4912339)
        gsm_id <- gsub("\\\\.telocal\\\\.cntTable\$", "", basename(f))
        
        # Read 2-column count table (feature ID and count)
        dt <- fread(f, header = FALSE)
        setnames(dt, c("feature", gsm_id))

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
    fwrite(merged_dt, file = "telocal_merged_counts.txt", sep = "\t", quote = FALSE)
    """
}

process EDGER_TELOCAL {
    tag "contrast_analysis_telocal"
    cpus 4
    memory '16 GB'
    publishDir "${params.output}/edger_telocal", mode: 'copy'

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

process SC_TE {
    tag { meta.gsm_id }
    cpus 5
    memory '80 GB'

    input:
    val meta

    output:
    path "${meta.gsm_id}_scTE.csv", emit: scte_dir

    script:
    """
    /work/InternalMedicine/s184335/sc//repos/scTE/bin/scTE \
        -i ${meta.bam} \
        -p ${task.cpus} \
        -x ${meta.scTE_idx} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${meta.gsm_id}_scTE
    """
}

process SC_TELOCAL {
    tag { meta.gsm_id }
    cpus 5
    memory '120 GB'

    input:
    val meta

    output:
    path "${meta.gsm_id}_scTEtx.csv", emit: scte_dir

    script:
    """
    /work/InternalMedicine/s184335/sc//repos/scTE/bin/scTE \
        -i ${meta.bam} \
        -p ${task.cpus} \
        -x ${meta.scTE_tx_idx} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${meta.gsm_id}_scTEtx
    """
}

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

process EDGER_SCTELOCAL {
    tag "contrast_analysis_scte"
    cpus 4
    memory '16 GB'
    publishDir "${params.output}/edger_scte", mode: 'copy'

    input:
    path counts_matrix
    path contrast_sheet

    output:
    path "*_scte_edger_results.csv"  , emit: results
    path "all_scte_edger_results.csv", emit: combined_results
    path "*_scte_mds.pdf"            , emit: mds_plots

    script:
    """
    #!/usr/bin/env Rscript

    library(edgeR)
    library(ggplot2)
    library(data.table)

    # 1. Load inputs
    meta_df <- read.table("${contrast_sheet}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    counts  <- read.table("${counts_matrix}", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

    # 2. Exclude ENSEMBL coding genes (ENSG for Human, ENSMUSG for Mouse)
    is_ensembl_gene <- grepl("^ENSG|^ENSMUSG", rownames(counts), ignore.case = TRUE)
    counts <- counts[!is_ensembl_gene, , drop = FALSE]

    if (nrow(counts) == 0) {
        stop("No scTE features remaining after removing ENSEMBL coding genes.")
    }

    ctrl_patterns <- c("wt", "ctrl", "control")
    reserved_cols <- c("gsm_id", "sample", "bam", "fastq", "fastq_1", "fastq_2", "tecount_file", "telocal_file", "scte_file")

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
        
        # Filter low-count features to optimize memory and runtime
        keep <- filterByExpr(dge)
        dge <- dge[keep, , keep.lib.sizes = FALSE]

        dge <- calcNormFactors(dge)

        # Output MDS Plot
        pdf(paste0(col_name, "_scte_mds.pdf"))
        plotMDS(dge, labels = sub_meta\$target_factor, main = paste("MDS (scTE):", col_name))
        dev.off()

        # 4. Fit quasi-likelihood GLM model
        design <- model.matrix(~ target_factor, data = sub_meta)
        dge <- estimateDisp(dge, design)
        fit <- glmQLFit(dge, design)

        # Extract contrast results for each treatment group
        for (i in seq_along(treatments)) {
            trt_group <- treatments[i]
            
            coef_idx <- paste0("target_factortrt_", trt_group)
            if (!coef_idx %in% colnames(design)) {
                coef_idx <- i + 1
            }

            qlf <- glmQLFTest(fit, coef = coef_idx)
            res_table <- topTags(qlf, n = Inf)\$table

            # Convert to data.table
            res_dt <- as.data.table(res_table, keep.rownames = "scte_id")

            # Standardize column names to log2FoldChange, pvalue, padj
            setnames(res_dt, 
                     old = c("logFC", "PValue", "FDR"), 
                     new = c("log2FoldChange", "pvalue", "padj"), 
                     skip_absent = TRUE)

            res_dt[, trt := trt_group]
            res_dt[, contrast_column := col_name]

            # Reorder columns with scte_id, trt, and contrast_column first
            setcolorder(res_dt, c("scte_id", "trt", "contrast_column"))

            output_name <- paste0(col_name, "_", trt_group, "_vs_ctrl_scte_edger_results.csv")
            fwrite(res_dt, file = output_name)

            all_results_list[[length(all_results_list) + 1]] <- res_dt
        }
    }

    # 5. Save combined table
    if (length(all_results_list) > 0) {
        combined_dt <- rbindlist(all_results_list, fill = TRUE)
        fwrite(combined_dt, file = "all_scte_edger_results.csv")
    }
    """
}

process IRFINDER_FASTQ {
    tag "${meta.gsm_id}"
    cpus 5
    memory '40 GB'

    input:
    val meta

    output:
    path "ir_out_${meta.gsm_id}",       emit: ir_dir
    path "ir_out_${meta.gsm_id}/*.txt", emit: ir_results

    script:
    def reads = meta.trim_r2 ? "${meta.trim_r1} ${meta.trim_r2}" : "${meta.trim_r1}"

    """
    IRFinder -m FASTQ \
        -r ${meta.irfinder_index} \
        -d ir_out_${meta.gsm_id} \
        -t ${task.cpus} \
        ${reads}
    """
}

process MERGE_IRFINDER {
    tag "merge_irfinder"
    cpus 2
    memory '16 GB'
    publishDir "${params.output}/irfinder_matrix", mode: 'copy'

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

process ALIGN_DNA {
    tag { meta.gsm_id }
    cpus 8
    memory 32.GB

    input:
    val meta 

    output:
    tuple val(meta), path("${meta.gsm_id}.sorted.bam"), path("${meta.gsm_id}.sorted.bam.bai")

    script:
    // Determine if we are in PE or SE mode for Bowtie2
    def read_input = meta.trim_r2 ? "-1 ${meta.trim_r1} -2 ${meta.trim_r2}" : "-U ${meta.trim_r1}"
    
    """
    bowtie2 --threads ${task.cpus} \
            --very-sensitive-local \
            --rg-id ${meta.gsm_id} \
            --rg "SM:${meta.gsm_id}" \
            --rg "PL:ILLUMINA" \
            -x ${meta.bowtie2_index} \
            ${read_input} \
            2> ${meta.gsm_id}_bowtie2.log | \
    samtools sort -@ ${task.cpus} -o ${meta.gsm_id}.sorted.bam -

    samtools index ${meta.gsm_id}.sorted.bam
    """
}

process MARK_DUPLICATES {
    tag { meta.gsm_id }
    cpus 2
    memory 10.GB

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.md.bam"), path("${meta.gsm_id}.md.bam.bai"), path("${meta.gsm_id}.metrics.txt")

    script:
    """
    # Use the environment variable provided by the module
    java -Xmx8G -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        I=${meta.bam} \
        O=${meta.gsm_id}.md.bam \
        M=${meta.gsm_id}.metrics.txt \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT

    # Ensure index naming is consistent
    if [ ! -f "${meta.gsm_id}.md.bam.bai" ]; then
        mv ${meta.gsm_id}.md.bai ${meta.gsm_id}.md.bam.bai
    fi
    """
}

process FILTER_BAM {
    tag { meta.gsm_id }
    cpus 2
    memory 4.GB

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.filtered.bam"), path("${meta.gsm_id}.filtered.bam.bai")

    script:
    """
    # -q 30: Quality score >= 30
    # -F 4: Exclude unmapped reads
    # -F 256: Exclude secondary alignments
    # -F 1024: Exclude PCR/optical duplicates
    # -b: output BAM format
    
    samtools view -b -q 30 -F 1804 ${meta.dedup_bam} > ${meta.gsm_id}.filtered.bam
    
    # Index the filtered BAM
    samtools index ${meta.gsm_id}.filtered.bam
    """
}

process GENERATE_BIGWIG {
    tag { meta.gsm_id }
    cpus 4
    memory 16.GB

    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}.bw")

    script:
    """
    bamCoverage \
        --bam ${meta.filt_bam} \
        --outFileName ${meta.gsm_id}.bw \
        --outFileFormat bigwig \
        --numberOfProcessors ${task.cpus} \
        --normalizeUsing CPM \
        --binSize 10 
    """
}

process MACS3_CALLPEAK_NoCONTROL {
    tag { meta.gsm_id }
    cpus 4
    memory 16.GB
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}_peaks.narrowPeak"), path("${meta.gsm_id}_summits.bed"), path("${meta.gsm_id}_peaks.xls")

    script:
    // Safely check for PE or PE_plus_R3 modes
    def format = (meta.mode.startsWith("PE")) ? "BAMPE" : "BAM"
    """
    macs3 callpeak \
        -t ${meta.filt_bam} \
        -f ${format} \
        -g ${params.genomes[params.genome].genomeSize} \
        -n ${meta.gsm_id} \
        -q 0.01 \
        --outdir .
    """
}

process SAMTOOLS_FLAGSTAT {
    tag { meta.gsm_id }
    cpus 1
    memory 2.GB

    input:
    tuple val(meta), path(bam)

    output:
    path "${meta.gsm_id}.flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${meta.gsm_id}.flagstat.txt
    """
}

process MULTIQC {
    cpus 10
    memory 40.GB
    
    input:
    path qc_inputs

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}

