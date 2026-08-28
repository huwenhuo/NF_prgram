process MERGE_TECOUNTS {
    tag "all_samples_tecount"
    cpus 2
    memory '8 GB'
    publishDir "${params.output}/tecount_analysis", mode: 'copy'

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
