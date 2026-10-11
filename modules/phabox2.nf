// =============================================================================
//  Module: PHABOX2
//  All-in-one phage characterisation using PhaBOX v2 end_to_end task
//  Runs: phamer, phagcn, phatyp, cherry, phavip
// =============================================================================

process PHABOX2 {
    tag "${sample_id}"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/07_phabox2/${sample_id}",
        mode: 'copy'
    )

    input:
    tuple val(sample_id), path(polished)

    output:
    tuple val(sample_id), path("${sample_id}_phabox2/"), emit: results
    path "${sample_id}_phabox2/${sample_id}_summary.tsv", emit: summary, optional: true

    script:
    def db          = params.phabox2_db      ?: error("ERROR: --phabox2_db is required")
    def min_len     = params.phabox2_min_len ?: 1000
    def skip_phamer = params.phabox2_skip_phamer ?: "Y"
    """
    phabox2 \\
        --task     end_to_end \\
        --dbdir    ${db} \\
        --outpth   ${sample_id}_phabox2 \\
        --contigs  ${polished} \\
        --threads  ${task.cpus} \\
        --len      ${min_len} \\
        --skip     ${skip_phamer} \\
        2>&1 | tee ${sample_id}_phabox2.log

    cp ${sample_id}_phabox2.log ${sample_id}_phabox2/

    [ -d ${sample_id}_phabox2 ] || \\
        { echo "ERROR: PhaBOX2 produced no output for ${sample_id}"; exit 1; }
    """

    stub:
    """
    mkdir -p ${sample_id}_phabox2
    echo -e "sample\tlifestyle\thost\ttaxonomy" > ${sample_id}_phabox2/${sample_id}_summary.tsv
    """
}
