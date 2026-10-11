// =============================================================================
//  Module: DNAAPLER
//  Reorient phage genome to start at the terminase large subunit (terL)
//  Runs after Medaka polishing, before annotation and characterisation.
//
//  Why: Flye assembles circular genomes starting at an arbitrary position
//  and may assemble on either strand. Standardising orientation ensures:
//    - Consistent gene coordinates across all isolates
//    - Correct pangenome graph alignment
//    - pharokka genome maps all start at the same gene
//    - No artificial "inversion" signal in pangenome visualisations
//
//  Uses MMseqs2 to find terL against the dnaapler database.
//  If terL is not found, the sequence is left unchanged (safe fallback).
// =============================================================================

process DNAAPLER {
    tag "${sample_id}"
    label 'process_medium'

    publishDir(
        path: "${params.outdir}/05_dnaapler",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log')) "logs/${filename}"
            else                          filename
        }
    )

    input:
    tuple val(sample_id), path(polished)

    output:
    tuple val(sample_id), path("${sample_id}_reoriented.fasta"), emit: reoriented
    path "${sample_id}_dnaapler.log",                             emit: log

    script:
    def threads = task.cpus
    """
    echo "=== Dnaapler reorientation: ${sample_id} ===" | tee    ${sample_id}_dnaapler.log
    echo "  Input  : ${polished}"                        | tee -a ${sample_id}_dnaapler.log
    echo "  Threads: ${threads}"                         | tee -a ${sample_id}_dnaapler.log
    echo ""                                               | tee -a ${sample_id}_dnaapler.log

    set +e
    dnaapler phage \\
        -i  ${polished} \\
        -o  dnaapler_out \\
        -p  ${sample_id} \\
        -t  ${threads} \\
        >> ${sample_id}_dnaapler.log 2>&1
    EXIT_CODE=\$?
    set -e

    cat ${sample_id}_dnaapler.log

    if [[ \$EXIT_CODE -eq 0 ]] && [[ -f "dnaapler_out/${sample_id}_reoriented.fasta" ]]; then
        cp dnaapler_out/${sample_id}_reoriented.fasta ${sample_id}_reoriented.fasta
        echo "" | tee -a ${sample_id}_dnaapler.log
        grep "Re-orienting\\|already starts\\|not found\\|coordinate\\|strand" \\
            ${sample_id}_dnaapler.log | tail -5 | tee -a ${sample_id}_dnaapler.log || true
    else
        # Fallback: if dnaapler fails, use the polished assembly as-is
        echo "WARNING: dnaapler failed — using polished assembly unchanged" \\
            | tee -a ${sample_id}_dnaapler.log
        cp ${polished} ${sample_id}_reoriented.fasta
    fi
    """

    stub:
    """
    cp ${polished} ${sample_id}_reoriented.fasta
    echo "stub" > ${sample_id}_dnaapler.log
    """
}
