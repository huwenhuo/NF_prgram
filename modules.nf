process DOWNLOAD_FASTQ {
    tag "${gsm_id}"
    
    input:
    val gsm_id

    output:
    tuple val(gsm_id), path("${gsm_id}*.fastq.gz", arity: '1..*')

    script:
    """
    # 1. Retrieve the SRR ID
    SRR_ID=\$(esearch -db sra -query "${gsm_id}" | efetch -format runinfo | grep "SRR" | cut -d',' -f1 | head -n 1)
    
    # 2. Prefetch the data (Robust download)
    prefetch \$SRR_ID
    
    # 3. Extract using the local folder (More reliable than web streaming)
    fasterq-dump --split-3 --include-technical --threads 4 \$SRR_ID/\$SRR_ID.sra
    
    # 4. Rename files
    for file in \$SRR_ID*; do
        suffix=\${file#\$SRR_ID}
        mv "\$file" "${gsm_id}\${suffix}"
    done
    
    # 5. Gzip the files
    gzip ${gsm_id}*.fastq
    """
}

process TRIM_FASTQ_PE {
    tag { gsm_id }
    cpus 4
    memory 8.GB

    input:
    tuple val(gsm_id), path(reads)

    output:
    tuple val(gsm_id), path("${gsm_id}_R1.trimmed.fastq.gz"), path("${gsm_id}_R2.trimmed.fastq.gz")

    script:
    """
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${gsm_id}_R1.trimmed.fastq.gz -O ${gsm_id}_R2.trimmed.fastq.gz \
        --disable_quality_filtering \
        --length_required 20 \
        --detect_adapter_for_pe \
        --thread ${task.cpus} \
        --html ${gsm_id}_fastp.html \
        --json ${gsm_id}_fastp.json
    """
}

process TRIM_FASTQ_SE {
    tag { gsm_id }
    cpus 4
    memory 8.GB

    input:
    tuple val(gsm_id), path(reads)

    output:
    tuple val(gsm_id), path("${gsm_id}_R1.trimmed.fastq.gz")

    script:
    """
    fastp \
        -i ${reads[0]} \
        -o ${gsm_id}_R1.trimmed.fastq.gz \
        --length_required 36 \
        --thread ${task.cpus} \
        --html ${gsm_id}_fastp.html \
        --json ${gsm_id}_fastp.json
    """
}

process ALIGN_DNA_PE {
    tag { gsm_id }
    cpus 8
    memory 32.GB

    input:
    tuple val(gsm_id), path(reads)
    val index_path

    output:
    tuple val(gsm_id), path("${gsm_id}.sorted.bam"), path("${gsm_id}.sorted.bam.bai")

    script:
    """
    bowtie2 --threads ${task.cpus} \
            --very-sensitive-local \
            -x ${index_path} \
            -1 ${reads[0]} -2 ${reads[1]} \
            2> ${gsm_id}_bowtie2.log | \
    samtools sort -@ ${task.cpus} -o ${gsm_id}.sorted.bam -

    samtools index ${gsm_id}.sorted.bam
    """
}

process ALIGN_DNA_SE {
    tag { gsm_id }
    cpus 8
    memory 32.GB

    input:
    tuple val(gsm_id), path(reads)
    val index_path

    output:
    tuple val(gsm_id), path("${gsm_id}.sorted.bam"), path("${gsm_id}.sorted.bam.bai")

    script:
    """
    bowtie2 --threads ${task.cpus} \
            --very-sensitive-local \
            -x ${index_path} \
            -U ${reads[0]} \
            2> ${gsm_id}_bowtie2.log | \
    samtools sort -@ ${task.cpus} -o ${gsm_id}.sorted.bam -

    samtools index ${gsm_id}.sorted.bam
    """
}

process ALIGN_RNA_STAR_PE {
    tag { gsm_id }
    cpus 16
    memory 64.GB

    input:
    tuple val(gsm_id), path(reads)
    val star_index

    output:
    tuple val(gsm_id), path("${gsm_id}.Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${star_index} \
         --readFilesIn ${reads[0]} ${reads[1]} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${gsm_id}. \
         --outStd Log
    """
}

process ALIGN_RNA_STAR_SE {
    tag { gsm_id }
    cpus 16
    memory 64.GB

    input:
    tuple val(gsm_id), path(reads)
    val star_index

    output:
    tuple val(gsm_id), path("${gsm_id}.Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${star_index} \
         --readFilesIn ${reads[0]} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${gsm_id}. \
         --outStd Log
    """
}

process MARK_DUPLICATES {
    tag { gsm_id }
    cpus 2
    memory 10.GB

    input:
    tuple val(gsm_id), path(sorted_bam), path(bai)

    output:
    tuple val(gsm_id), path("${gsm_id}.md.bam"), path("${gsm_id}.md.bam.bai"), path("${gsm_id}.metrics.txt")

    script:
    """
    # Use the environment variable provided by the module
    java -Xmx8G -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        I=${sorted_bam} \
        O=${gsm_id}.md.bam \
        M=${gsm_id}.metrics.txt \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT

    # Ensure index naming is consistent
    if [ ! -f "${gsm_id}.md.bam.bai" ]; then
        mv ${gsm_id}.md.bai ${gsm_id}.md.bam.bai
    fi
    """
}

process FILTER_BAM {
    tag { gsm_id }
    cpus 2
    memory 4.GB

    input:
    tuple val(gsm_id), path(bam), path(bai), path(metrics)

    output:
    tuple val(gsm_id), path("${gsm_id}.filtered.bam"), path("${gsm_id}.filtered.bam.bai")

    script:
    """
    # -q 30: Quality score >= 30
    # -F 4: Exclude unmapped reads
    # -F 256: Exclude secondary alignments
    # -F 1024: Exclude PCR/optical duplicates
    # -b: output BAM format
    
    samtools view -b -q 30 -F 1804 ${bam} > ${gsm_id}.filtered.bam
    
    # Index the filtered BAM
    samtools index ${gsm_id}.filtered.bam
    """
}

process MACS3_CALLPEAK {
    tag { gsm_id }
    cpus 4
    memory 16.GB

    input:
    // We expect the treatment and control BAMs to be paired
    tuple val(gsm_id), path(treatment_bam), path(treatment_bai), path(control_bam), path(control_bai)

    output:
    tuple val(gsm_id), path("${gsm_id}_peaks.narrowPeak"), path("${gsm_id}_summits.bed"), path("${gsm_id}_peaks.xls")

    script:
    // -f BAMPE is for paired-end; use -f BAM if single-end
    // -g hs is for Human (GRCh38), mm for Mouse (mm10)
    // -q 0.01 is the default q-value cutoff
    """
    macs3 callpeak \
        -t ${treatment_bam} \
        -c ${control_bam} \
        -f BAMPE \
        -g ${params.genomes[params.genome].genomeSize} \
        -n ${gsm_id} \
        -q 0.01 \
        --outdir .
    """
}

process MACS3_CALLPEAK_NoCONTROL {
    tag { gsm_id }
    cpus 4
    memory 16.GB

    input:
    tuple val(gsm_id), path(bam), path(bai)

    output:
    tuple val(gsm_id), path("${gsm_id}_peaks.narrowPeak"), path("${gsm_id}_summits.bed"), path("${gsm_id}_peaks.xls")

    script:
    """
    macs3 callpeak \
        -t ${bam} \
        -f BAMPE \
        -g ${params.genomes[params.genome].genomeSize} \
        -n ${gsm_id} \
        -q 0.01 \
        --outdir .
    """
}

process GENERATE_BIGWIG {
    tag { gsm_id }
    cpus 4
    memory 16.GB

    input:
    tuple val(gsm_id), path(bam), path(bai)

    output:
    tuple val(gsm_id), path("${gsm_id}.bw")

    script:
    """
    bamCoverage \
        --bam ${bam} \
        --outFileName ${gsm_id}.bw \
        --outFileFormat bigwig \
        --numberOfProcessors ${task.cpus} \
        --normalizeUsing CPM \
        --binSize 10 \
        --extendReads
    """
}


