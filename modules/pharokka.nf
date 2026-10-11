// =============================================================================
//  Module: PHAROKKA
//  Comprehensive phage genome annotation
//
//  Input  : polished FASTA from MEDAKA (same input as PHABOX2/ASSEMBLY_STATS)
//  DB     : params.pharokka_db — path to pharokka databases
//  Outputs: GenBank (.gbk), GFF (.gff), protein FASTA (.faa),
//           nucleotide genes (.ffn), summary TSV
//
//  Gene predictor: phanotate (recommended for phage genomes)
// =============================================================================

process PHAROKKA {
    tag "${sample_id}"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/08_pharokka/${sample_id}",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log')) "logs/${filename}"
            else                          filename
        }
    )

    input:
    tuple val(sample_id), path(polished)

    output:
    tuple val(sample_id), path("pharokka_out/"),                          emit: results
    tuple val(sample_id), path("pharokka_out/${sample_id}.gbk"),          emit: genbank,  optional: true
    tuple val(sample_id), path("pharokka_out/${sample_id}.gff"),          emit: gff,      optional: true
    tuple val(sample_id), path("pharokka_out/${sample_id}.faa"),          emit: proteins, optional: true
    path "${sample_id}_pharokka.log",                                      emit: log

    script:
    def db      = params.pharokka_db ?: error("ERROR: --pharokka_db is required")
    def threads = task.cpus
    """
    # Print version for diagnostics
    pharokka.py -V >> ${sample_id}_pharokka.log 2>&1 || true

    echo "=== Pharokka annotation: ${sample_id} ===" | tee -a ${sample_id}_pharokka.log
    echo "  Input  : ${polished}"                    | tee -a ${sample_id}_pharokka.log
    echo "  DB     : ${db}"                          | tee -a ${sample_id}_pharokka.log
    echo "  Threads: ${threads}"                     | tee -a ${sample_id}_pharokka.log
    echo ""                                           | tee -a ${sample_id}_pharokka.log

    # PYTHONWARNINGS=ignore prevents phanotate pkg_resources deprecation warning
    # from being written into GFF output and corrupting it.
    # --mmseqs2_only skips PyHMMER which crashes due to API changes in newer pyhmmer.
    set +e
    PYTHONWARNINGS=ignore pharokka.py \\
        -i  ${polished} \\
        -o  pharokka_out \\
        -d  ${db} \\
        -t  ${threads} \\
        -g  phanotate \\
        -p  ${sample_id} \\
        -f  \\
        --mmseqs2_only \\
        >> ${sample_id}_pharokka.log 2>&1
    EXIT_CODE=\$?
    set -e

    cat ${sample_id}_pharokka.log

    [ \$EXIT_CODE -eq 0 ] || \\
        { echo "ERROR: Pharokka failed for ${sample_id} (exit \$EXIT_CODE)"; exit 1; }

    [ -d pharokka_out ] || \\
        { echo "ERROR: Pharokka produced no output for ${sample_id}"; exit 1; }

    cp ${sample_id}_pharokka.log pharokka_out/
    ls pharokka_out/ >> ${sample_id}_pharokka.log
    """

    stub:
    """
    mkdir -p pharokka_out
    echo "stub gbk"  > pharokka_out/${sample_id}.gbk
    echo "stub gff"  > pharokka_out/${sample_id}.gff
    echo "stub faa"  > pharokka_out/${sample_id}.faa
    echo "stub"      > ${sample_id}_pharokka.log
    """
}
