// =============================================================================
//  Module: MLST
//  Sequence typing of the host assembly (tseemann/mlst).
//  params.mlst_scheme = null → automatic scheme detection.
// =============================================================================

process MLST {
    tag "${sample_id}"
    label 'process_low'

    publishDir(
        path: "${params.outdir}/09_host/${sample_id}",
        mode: 'copy'
    )

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_mlst.tsv"), emit: results

    script:
    def scheme = params.mlst_scheme ? "--scheme ${params.mlst_scheme}" : ''
    """
    mlst ${scheme} ${assembly} > ${sample_id}_mlst.tsv
    """

    stub:
    """
    printf '${sample_id}_host_polished.fasta\\tklebsiella\\t3691\\tgapA(2)\\tinfB(1)\\tmdh(1)\\tpgi(1)\\tphoE(1)\\trpoB(1)\\ttonB(1)\\n' > ${sample_id}_mlst.tsv
    """
}
