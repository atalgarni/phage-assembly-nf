// =============================================================================
//  Module: PORECHOP_ABI
//  Adapter detection and trimming for ONT reads
// =============================================================================

process PORECHOP_ABI {
    tag "${sample_id}"
    label 'process_medium'

    publishDir(
        path: "${params.outdir}/01_porechop_abi",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log'))     "logs/${filename}"
            else                              filename
        }
    )

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"),  emit: trimmed_reads
    path "${sample_id}_porechop.log",                              emit: log

    script:
    def abi_flag       = params.abi            ? '-abi'                  : ''
    def no_split_flag  = params.no_split       ? '--no_split'            : ''
    def discard_mid    = params.discard_middle ? '--discard_middle'      : ''
    def min_split      = "--min_split_read_size ${params.min_split_read_size}"

    """
    porechop_abi \\
        ${abi_flag} \\
        -i  ${reads} \\
        -o  ${sample_id}_trimmed.fastq.gz \\
        --format fastq.gz \\
        -t  ${task.cpus} \\
        ${no_split_flag} \\
        ${discard_mid} \\
        ${min_split} \\
        2>&1 | tee ${sample_id}_porechop.log

    # Fail loudly if the output is empty
    [ -s ${sample_id}_trimmed.fastq.gz ] || \\
        { echo "ERROR: Porechop_ABI produced no output for ${sample_id}"; exit 1; }
    """

    stub:
    // Stub for fast pipeline testing without running the real tool
    """
    echo "Stub: Porechop_ABI for ${sample_id}"
    gzip -c /dev/null > ${sample_id}_trimmed.fastq.gz
    echo "stub log" > ${sample_id}_porechop.log
    """
}
