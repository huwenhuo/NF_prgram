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
    tag { meta.gsm_id }
    cpus 10
    memory '60 GB'
    
    input:
    val meta

    output:
    tuple val(meta), path("${meta.gsm_id}_TE_Aligned.sortedByCoord.out.bam"), path("${meta.gsm_id}_TE_Aligned.sortedByCoord.out.bam.bai")

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

    samtools index ${meta.gsm_id}_TE_Aligned.sortedByCoord.out.bam
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

    # Resolve absolute paths for the reference directory
    INDEX_ABS=\$(readlink -f "${meta.irfinder_index}")

    # Resolve absolute path(s) for input trimmed reads
    READ_FILES=""
    for f in ${reads}; do
        READ_FILES="\${READ_FILES} \$(readlink -f \$f)"
    done

    singularity exec --bind /project,/archive,/home,/endosome,/work,\$PWD ${meta.img_irfinder} \
        IRFinder -m FASTQ \
        -r \${INDEX_ABS} \
        -d ir_out_${meta.gsm_id} \
        -t ${task.cpus} \
        \${READ_FILES}
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

