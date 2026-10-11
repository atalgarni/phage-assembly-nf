// =============================================================================
//  Module: HOST_MEDAKA
//  Polish the host assembly; falls back to the Flye draft if Medaka fails.
// =============================================================================

process HOST_MEDAKA {
    tag "${sample_id}"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/09_host/${sample_id}",
        mode: 'copy',
        saveAs: { filename -> filename.endsWith('.log') ? "logs/${filename}" : filename }
    )

    input:
    tuple val(sample_id), path(host_reads), path(draft)

    output:
    tuple val(sample_id), path("${sample_id}_host_polished.fasta"), emit: polished
    path "${sample_id}_host_medaka.log",                             emit: log

    script:
    def model = params.medaka_model ?: "r1041_e82_400bps_hac_v4.3.0"
    """
    set +e
    medaka_consensus \\
        -i ${host_reads} \\
        -d ${draft} \\
        -o medaka_host \\
        -m ${model} \\
        -t ${task.cpus} \\
        > ${sample_id}_host_medaka.log 2>&1
    EXIT_CODE=\$?
    set -e

    if [[ \$EXIT_CODE -eq 0 ]] && [[ -s medaka_host/consensus.fasta ]]; then
        cp medaka_host/consensus.fasta ${sample_id}_host_polished.fasta
    else
        echo "WARNING: Medaka failed for host of ${sample_id} — using Flye draft" \\
            | tee -a ${sample_id}_host_medaka.log
        cp ${draft} ${sample_id}_host_polished.fasta
    fi
    """

    stub:
    """
    cp ${draft} ${sample_id}_host_polished.fasta
    echo stub > ${sample_id}_host_medaka.log
    """
}
