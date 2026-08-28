process PATHWAY_ANALYSIS {
    tag "fgsea_${contrast_file.baseName}"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/coding_gene", mode: 'copy'

    input:
    path contrast_file

    output:
    path "*_fgsea_hallmark.csv"   , emit: hallmark_results, optional: true
    path "*_fgsea_gobp.csv"       , emit: gobp_results    , optional: true
    path "*_fgsea_summary.pdf"    , emit: plots           , optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(fgsea)
    library(msigdbr)
    library(ggplot2)
    library(ggpubr)
    library(AnnotationDbi)

    # 1. Read DESeq2 result CSV
    res_dt <- fread("${contrast_file}")

    if (!"gene_id" %in% names(res_dt) || !"stat" %in% names(res_dt)) {
        stop("Input file missing required 'gene_id' or 'stat' columns.")
    }

    # Clean ENSEMBL IDs
    res_dt[, clean_gene_id := sub("\\\\..*", "", gene_id)]

    # Remove NA or infinite test statistics
    res_dt <- res_dt[!is.na(stat) & is.finite(stat)]

    if (nrow(res_dt) < 100) {
        message("Skipping GSEA: Too few valid genes in input table.")
        q(save = "no", status = 0)
    }

    # Detect species automatically from ENSEMBL ID prefix
    sample_gene <- res_dt\$clean_gene_id[1]
    if (grepl("^ENSG", sample_gene, ignore.case = TRUE)) {
        species_name <- "Homo sapiens"
        org_db_pkg   <- "org.Hs.eg.db"
    } else if (grepl("^ENSMUSG", sample_gene, ignore.case = TRUE)) {
        species_name <- "Mus musculus"
        org_db_pkg   <- "org.Mm.eg.db"
    } else {
        stop("Could not infer species from ENSEMBL IDs. Expected ENSG or ENSMUSG.")
    }

    library(org_db_pkg, character.only = TRUE)
    org_db <- get(org_db_pkg)

    # 2. Map ENSEMBL IDs to Entrez IDs
    entrez_map <- mapIds(
        org_db,
        keys      = res_dt\$clean_gene_id,
        column    = "ENTREZID",
        keytype   = "ENSEMBL",
        multiVals = "first"
    )

    res_dt[, entrez_id := entrez_map[clean_gene_id]]
    res_dt <- res_dt[!is.na(entrez_id) & !duplicated(entrez_id)]

    # 3. Create ranked gene vector using DESeq2 Wald statistic
    ranks <- res_dt\$stat
    names(ranks) <- res_dt\$entrez_id
    ranks <- sort(ranks, decreasing = TRUE)

    # 4. Helper function to fetch gene sets and run fgsea
    run_fgsea_collection <- function(category, subcategory = NULL) {
        msig_df <- msigdbr(species = species_name, category = category, subcategory = subcategory)
        pathways <- split(x = msig_df\$entrez_gene, f = msig_df\$gs_name)
        
        fgsea_res <- fgsea(
            pathways = pathways,
            stats    = ranks,
            minSize  = 15,
            maxSize  = 500,
            nperm    = 10000
        )
        return(as.data.table(fgsea_res))
    }

    base_prefix <- sub("_coding_deseq2_results.*", "", "${contrast_file}")

    # 5. Run GSEA on Hallmark sets (H)
    hallmark_res <- run_fgsea_collection(category = "H")
    if (nrow(hallmark_res) > 0) {
        hallmark_res <- hallmark_res[order(padj, pval)]
        hallmark_res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
        fwrite(hallmark_res, file = paste0(base_prefix, "_fgsea_hallmark.csv"))

        # --- DEBUG: Print p53 Pathway status in process log ---
        p53_row <- hallmark_res[grepl("P53_PATHWAY", pathway, ignore.case = TRUE)]
        if (nrow(p53_row) > 0) {
            p53_rank <- which(hallmark_res\$pathway == p53_row\$pathway[1])
            message(sprintf("=== [DEBUG] P53 Pathway Status in %s ===", base_prefix))
            message(sprintf("Rank: %d / %d | NES: %.3f | pval: %.4e | padj: %.4e", 
                            p53_rank, nrow(hallmark_res), p53_row\$NES[1], p53_row\$pval[1], p53_row\$padj[1]))
        } else {
            message("=== [DEBUG] P53 Pathway not found in Hallmark results ===")
        }
    }

    # 6. Run GSEA on GO Biological Processes (C5:BP)
    gobp_res <- run_fgsea_collection(category = "C5", subcategory = "BP")
    if (nrow(gobp_res) > 0) {
        gobp_res <- gobp_res[order(padj, pval)]
        gobp_res[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
        fwrite(gobp_res, file = paste0(base_prefix, "_fgsea_gobp.csv"))
    }

    # 7. Generate Dotplot AND Barplot for the exact same selected pathways
    if (nrow(hallmark_res) > 0 && sum(!is.na(hallmark_res\$padj)) > 0) {
        
        hallmark_res[, pathway_clean := gsub("^HALLMARK_", "", pathway)]

        # Priority pathways to retain
        priority_pathways <- c("P53_PATHWAY", "APOPTOSIS", "DNA_REPAIR", "MYC_TARGETS_V1")
        
        # Select top 20 overall pathways by significance
        top_ranked <- head(hallmark_res[order(padj, pval)], 20)\$pathway_clean

        # Unified pathway selection list
        pathway_sel <- unique(c(priority_pathways, top_ranked))
        
        plotdat <- hallmark_res[pathway_clean %in% pathway_sel, ]

        if (nrow(plotdat) > 0) {
            # Sort by NES for clean display
            plotdat <- plotdat[order(NES), ]
            plotdat[, pathway_clean := factor(pathway_clean, levels = plotdat\$pathway_clean)]
            plotdat[, sig := fifelse(padj < 0.1 & NES > 0, "Up", fifelse(padj < 0.1 & NES < 0, "Down", "NSig"))]
            
            # --- Plot 1: Original Dotplot ---
            p_dot <- ggplot(plotdat, aes(x = NES, y = pathway_clean)) +
                geom_point(aes(size = -log10(padj), color = NES)) +
                scale_color_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
                theme_minimal(base_size = 8) +
                labs(
                    title = paste("Hallmark Pathways (Dotplot):", base_prefix),
                    x = "Normalized Enrichment Score (NES)",
                    y = "",
                    size = "-log10(adj p-value)"
                ) +
                theme(
                    axis.text.y = element_text(size = 7),
                    axis.text.x = element_text(size = 7),
                    plot.title  = element_text(hjust = 0.5)
                )

            # --- Plot 2: Directional Barplot ---
            sig_colors <- c("Up" = "#E41A1C", "NSig" = "#999999", "Down" = "#377EB8")

            p_bar <- ggplot(plotdat, aes(x = pathway_clean, y = NES, fill = sig)) +
                geom_col() +
                coord_flip() +
                scale_fill_manual(values = sig_colors) +
                labs(
                    x = "",
                    y = "Normalized Enrichment Score (NES)",
                    title = paste("Hallmark Pathways (Barplot):", base_prefix),
                    fill = "Significance (padj < 0.1)"
                ) +
                theme_pubr(base_size = 8) +
                theme(
                    axis.text.y = element_text(size = 7),
                    axis.text.x = element_text(size = 7),
                    plot.title  = element_text(hjust = 0.5)
                )

            # Save both pages to the PDF
            pdf_height <- max(4, nrow(plotdat) * 0.22)
            pdf(paste0(base_prefix, "_fgsea_summary.pdf"), width = 6.5, height = pdf_height)
            print(p_dot)
            print(p_bar)
            dev.off()
        }
    }
    """
}