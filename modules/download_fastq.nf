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
