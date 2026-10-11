// =============================================================================
//  Module: MEDAKA
//  Polish Flye draft assembly using medaka_consensus wrapper
//  Model: r1041_e82_400bps_hac_v4.3.0
//         (R10.4.1 / e8.2 pore / 400bps / Dorado HAC v4.3.0)
// =============================================================================

process MEDAKA {
    tag "${sample_id}"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/04_medaka/${sample_id}",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log')) "logs/${filename}"
            else                          filename
        }
    )

    input:
    tuple val(sample_id), path(reads), path(draft)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), emit: polished
    path "${sample_id}_medaka.log",                            emit: log

    script:
    def model   = params.medaka_model ?: "r1041_e82_400bps_hac_v4.3.0"
    def threads = task.cpus
    """
    echo "=== Medaka polishing: ${sample_id} ===" | tee    ${sample_id}_medaka.log
    echo "  Draft  : ${draft}"                    | tee -a ${sample_id}_medaka.log
    echo "  Reads  : ${reads}"                    | tee -a ${sample_id}_medaka.log
    echo "  Model  : ${model}"                    | tee -a ${sample_id}_medaka.log
    echo "  Threads: ${threads}"                  | tee -a ${sample_id}_medaka.log
    echo ""                                        | tee -a ${sample_id}_medaka.log

    medaka_consensus \\
        -i ${reads} \\
        -d ${draft} \\
        -o medaka_out \\
        -m ${model} \\
        -t ${threads} \\
        2>&1 | tee -a ${sample_id}_medaka.log

    cp medaka_out/consensus.fasta ${sample_id}_polished.fasta

    [ -s ${sample_id}_polished.fasta ] || \\
        { echo "ERROR: Medaka produced no output for ${sample_id}"; exit 1; }

    echo "" | tee -a ${sample_id}_medaka.log
    echo "Polished contigs:" | tee -a ${sample_id}_medaka.log
    grep "^>" ${sample_id}_polished.fasta | tee -a ${sample_id}_medaka.log
    """

    stub:
    """
    echo ">contig_1_polished" > ${sample_id}_polished.fasta
    echo "ATCGATCGATCG"      >> ${sample_id}_polished.fasta
    echo "stub"               > ${sample_id}_medaka.log
    """
}
