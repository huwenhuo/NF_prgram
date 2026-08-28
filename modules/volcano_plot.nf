process VOLCANO_PLOT {
    tag "volcano_${contrast_name}"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/coding_gene", mode: 'copy'

    input:
    tuple val(contrast_name), path(results_csv)

    output:
    path "*_volcano.pdf", emit: pdf, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(ggplot2)
    library(ggrepel)

    res_dt <- fread("${results_csv}")

    id_col <- intersect(c("symbol", "gene_id", "tecount"), colnames(res_dt))[1]
    if (is.na(id_col)) id_col <- colnames(res_dt)[1]

    res_dt[, sig := "NS"]
    res_dt[padj < 0.05 & log2FoldChange > 1, sig := "Upregulated"]
    res_dt[padj < 0.05 & log2FoldChange < -1, sig := "Downregulated"]
    res_dt[, sig := factor(sig, levels = c("Upregulated", "Downregulated", "NS"))]

    top_labels <- head(res_dt[sig != "NS"][order(padj)], 10)

    cols <- c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "NS" = "#999999")

    p <- ggplot(res_dt, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_manual(values = cols) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
        labs(
            title = paste("Volcano Plot:", "${contrast_name}"),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted P-Value",
            color = "Significance"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right"
        )

    if (nrow(top_labels) > 0) {
        p <- p + geom_text_repel(
            data = top_labels,
            aes(label = get(id_col)),
            size = 3,
            box.padding = 0.3,
            max.overlaps = 20,
            show.legend = FALSE
        )
    }

    pdf("${contrast_name}_volcano.pdf", width = 6, height = 5)
    print(p)
    dev.off()
    """
}