nextflow.enable.dsl=2

// ---- Parameters ----
params.samplesheet = "samplesheet.csv"
params.fastq_dir = "./"
params.outdir = "results"

def genome_params = params.genomes[params.genome]

println "STAR index: ${genome_params.star_index}"
println "FASTA: ${genome_params.fasta}"
println "GTF: ${genome_params.gtf}"


process DOWNLOAD_FASTQ {
    tag "$run_id"

    input:
    val run_id

    output:
    tuple val(run_id), path("${run_id}*.fastq")

    script:
    """
    module load sratoolkit/3.0.0
    fasterq-dump ${run_id} --split-files
    """
}

process TRIM_FASTQ {
    tag "${sample_id}"

    cpus = 2
    memory = 10.GB
    time = 10.h

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*.trimmed.fastq.gz")

    script:
    if (reads instanceof List && reads.size() == 2) {
        """
        module load singularityce/4.1.0
        
        R1_ABS=\$(readlink -f "${reads[0]}")
        R2_ABS=\$(readlink -f "${reads[1]}")

        singularity exec --bind /project,/archive,/home,/endosome,/work,\$PWD ${params.img_fastp} fastp \
            -i \Professional \${R1_ABS} \
            -I \${R2_ABS} \
            -o ${sample_id}_R1.trimmed.fastq.gz \
            -O ${sample_id}_R2.trimmed.fastq.gz \
            --detect_adapter_for_pe \
            --length_required 36 \
            --thread ${task.cpus} \
            --html ${sample_id}_fastp.html \
            --json ${sample_id}_fastp.json
        """
    } else {
        """
        module load singularityce/4.1.0
        
        SINGLE_READ=\$(echo "${reads}" | awk '{print \$1}')
        R1_ABS=\$(readlink -f "\${SINGLE_READ}")

        singularity exec --bind /project,/archive,/home,/endosome,/work,\$PWD ${params.img_fastp} fastp \
            -i \${R1_ABS} \
            -o ${sample_id}_R1.trimmed.fastq.gz \
            --length_required 36 \
            --thread ${task.cpus} \
            --html ${sample_id}_fastp.html \
            --json ${sample_id}_fastp.json
        """
    }
}


process STAR_teALIGNMENT {
    tag "${sample_id}"

    cpus = 10
    memory = 60.GB

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_te_Aligned.sortedByCoord.out.bam"), path("${sample_id}_te_Aligned.sortedByCoord.out.bam.bai")

    script:
    def read_files = (reads instanceof List) ? reads.join(" ") : reads
    """
    module load star/2.7.11b
    STAR --genomeDir ${genome_params.star_index}  \
         --runThreadN ${task.cpus} \
         --runMode alignReads \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 1000 \
         --outSAMmultNmax -1 \
         --outMultimapperOrder Random \
         --winAnchorMultimapNmax 1000 \
         --alignTranscriptsPerReadNmax 1000 \
         --alignMatesGapMax 350 \
         --readFilesIn ${read_files} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_te_ 
    module load samtools/1.22.1
    samtools index ${sample_id}_te_Aligned.sortedByCoord.out.bam
    """
}

process STAR_ALIGNMENT {
    tag "${sample_id}"

    publishDir "${params.outdir}/gene_counts", mode: 'copy', pattern: "*_ReadsPerGene.out.tab"

    cpus = 8
    memory = 60.GB

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_ReadsPerGene.out.tab"), path("${sample_id}_Aligned.sortedByCoord.out.bam")

    script:
    def read_files = (reads instanceof List) ? reads.join(" ") : reads
    """
    module load star/2.7.11b
    module load samtools/1.22.1
    STAR --genomeDir ${genome_params.star_index}  \
         --runThreadN ${task.cpus} \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 10 \
         --quantMode GeneCounts \
         --readFilesIn ${read_files} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}_ 
    samtools index ${sample_id}_Aligned.sortedByCoord.out.bam
    """
}


process TECOUNT {
    tag "${sample_id}"
    publishDir "${params.outdir}/te_count", mode: 'copy'

    cpus = 2
    memory = 20.GB
    time = 20.h

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(GTF_FILE), path(TE_GTF_FILE)

    output:
    tuple val(sample_id), path("${sample_id}.tecount.cntTable")

    script:
    """
    module load singularityce/4.1.0
    singularity exec ${params.img_tecount} TEcount \
        --sortByPos --format BAM --mode multi \
        -b ${bam_file} \
        --GTF ${GTF_FILE} \
        --TE ${TE_GTF_FILE} \
        --project ${sample_id}.tecount
    touch ${sample_id}.tecount.cntTable
    """
}

process TELOCAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/telocal_count/", mode: 'copy'

    cpus = 1
    memory = 10.GB
    time = 20.h

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), path(GTF_FILE), path(TE_loc_GTF_FILE)

    output:
    tuple val(sample_id), path("${sample_id}.telocal.cntTable")

    script:
    """
    module load singularityce/4.1.0
    singularity exec ${params.img_telocal} TElocal \
        --sortByPos -b ${bam_file} \
        --GTF ${GTF_FILE} \
        --TE ${TE_loc_GTF_FILE} \
        --stranded reverse \
        --project ${sample_id}.telocal
    touch ${sample_id}.telocal.cntTable 
    """
}

process SC_TELOCAL {
    tag "${sample_id}"
    publishDir "${params.outdir}/scte_count/", mode: 'copy'

    cpus = 2
    memory = 30.GB
    maxForks = 2   

    input:
    tuple val(sample_id), path(bam), val(idx_file) 

    output:
    path "${sample_id}_scTEtx.csv", emit: scTE_dir

    script:
    """
    scTE \
        -i ${bam} \
        -p ${task.cpus} \
        -x ${idx_file} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${sample_id}_scTEtx
    """
}

process SC_TE {
    tag "${sample_id}"
    publishDir "${params.outdir}/scte_count/", mode: 'copy'

    cpus = 5
    memory = 40.GB

    input:
    tuple val(sample_id), path(bam), path(idx_file)

    output:
    path "${sample_id}_scTE.csv", emit: scTE_dir

    script:
    """
    scTE \
        -i ${bam} \
        -p ${task.cpus} \
        -x ${idx_file} \
        --hdf5 False \
        -CB False \
        -UMI False \
        -o ${sample_id}_scTE
    """
}


process IRFINDER_FASTQ {
    tag "${sample_id}"
    publishDir "${params.outdir}/irfinder", mode: 'copy'

    cpus = 5
    memory = 40.GB

    input:
    tuple val(sample_id), path(reads)
    path irfinder_index
    path img_irfinder

    output:
    path "ir_out_${sample_id}", emit: ir_dir
    path "ir_out_${sample_id}/*.txt", emit: ir_results

    script:
    """
    module load singularityce/4.1.0

    # Resolve absolute paths for the reference directory
    INDEX_ABS=\$(readlink -f "${irfinder_index}")

    # Handle fastq inputs cleanly for both single and paired end
    if [ \$(echo "${reads}" | wc -w) -eq 2 ]; then
        # Paired-end files
        R1_ABS=\$(readlink -f \$(echo "${reads}" | awk '{print \$1}'))
        R2_ABS=\$(readlink -f \$(echo "${reads}" | awk '{print \$2}'))
        READ_FILES="\${R1_ABS} \${R2_ABS}"
    else
        # Single-end file
        R1_ABS=\$(readlink -f \$(echo "${reads}" | awk '{print \$1}'))
        READ_FILES="\${R1_ABS}"
    fi

    singularity exec --bind /project,/archive,/home,/endosome,/work,\$PWD ${img_irfinder} \
        IRFinder -m FASTQ \
        -r \${INDEX_ABS} \
        -d ir_out_${sample_id} \
        -t ${task.cpus} \
        \${READ_FILES}
    """
}

workflow {

    // ----------------------
    // Step 0: Download SRA FastQs
    // ----------------------
    channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> row.Run }
        .buffer(size: 1)
	.flatten()
        .set { run_ids_ch }

    sra_fastq_ch = DOWNLOAD_FASTQ(run_ids_ch)

    // ----------------------
    // Step 1: Trim FASTQs (Will handle single or paired dynamically)
    // ----------------------
    fastq_ch = sra_fastq_ch | TRIM_FASTQ

    // ----------------------
    // Step 2a: STAR alignment for normal genes
    // ----------------------
    star_gene_bam_ch = fastq_ch | STAR_ALIGNMENT

    // ----------------------
    // Step 2b: STAR alignment for TE counting
    // ----------------------
    star_te_bam_ch = fastq_ch | STAR_teALIGNMENT

    // ----------------------
    // Step 5: TE SC
    // ----------------------
    star_te_bam_ch
        .map { sample_id, bam, bai ->
            tuple(sample_id, bam, file(genome_params.scTE_idx))
        }
        | SC_TE

    // ----------------------
    // Step 6: TE SC with transcript idx
    // ----------------------
    star_te_bam_ch
        .map { sample_id, bam, bai -> 
            tuple(sample_id, bam, genome_params.scTE_tx_idx)
        }
        | SC_TELOCAL

    // ----------------------
    // Step 7: Intron Retention Analysis via FASTQ
    // ----------------------
    IRFINDER_FASTQ(
        fastq_ch,
        genome_params.irfinder_index,
        params.img_irfinder
    )
}
