// =============================================================================
//  Module: KAPTIVE
//  Capsule (K) and O-antigen (O) locus typing of the host assembly (Kaptive 3).
//  Databases are Kaptive keywords, e.g. kpsc_k / kpsc_o (Klebsiella),
//  ab_k / ab_o (Acinetobacter baumannii).
// =============================================================================

process KAPTIVE {
    tag "${sample_id}"
    label 'process_medium'

    publishDir(
        path: "${params.outdir}/09_host/${sample_id}/kaptive",
        mode: 'copy'
    )

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_kaptive_K.tsv"), path("${sample_id}_kaptive_O.tsv"), emit: results
    path "*.png", emit: plots, optional: true

    script:
    """
    kaptive assembly ${params.kaptive_k_db} ${assembly} \\
        -o ${sample_id}_kaptive_K.tsv \\
        -p ./ --plot-fmt png \\
        -t ${task.cpus}

    kaptive assembly ${params.kaptive_o_db} ${assembly} \\
        -o ${sample_id}_kaptive_O.tsv \\
        -t ${task.cpus}
    """

    stub:
    """
    printf 'Assembly\\tBest match locus\\tBest match type\\tMatch confidence\\tProblems\\tIdentity\\n${sample_id}\\tKL46\\tK46\\tTypeable\\t\\t99.5%%\\n' > ${sample_id}_kaptive_K.tsv
    printf 'Assembly\\tBest match locus\\tBest match type\\tMatch confidence\\tProblems\\tIdentity\\n${sample_id}\\tOL3\\tO3\\tTypeable\\t\\t99.9%%\\n' > ${sample_id}_kaptive_O.tsv
    """
}
