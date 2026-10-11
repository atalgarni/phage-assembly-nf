// =============================================================================
//  Module: FILTLONG
//  Length-filter ONT reads before assembly.
//  --min_length  : discard reads shorter than this (default 1000 bp)
//  --target_bases: subsample to this many bases, longest-first (optional)
//                  useful for very high coverage samples (e.g. host contamination)
//                  e.g. 11000000 = 100x coverage of a 110 kb phage
// =============================================================================

process FILTLONG {
    tag "${sample_id}"
    label 'process_low'

    publishDir(
        path: "${params.outdir}/02_filtlong",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log')) "logs/${filename}"
            else                          filename
        }
    )

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered_reads
    path "${sample_id}_filtlong.log",                              emit: log

    script:
    def min_len      = params.min_length   ?: 1000
    // Only add --target_bases if explicitly set (not null/empty)
    def target_arg   = params.target_bases ? "--target_bases ${params.target_bases}" : ''
    """
    # gunzip -c used instead of zcat: macOS zcat requires .Z extension
    BEFORE=\$(gunzip -c ${reads} | grep -c '^@' || true)

    filtlong \\
        --min_length ${min_len} \\
        ${target_arg} \\
        ${reads} \\
        2> ${sample_id}_filtlong.log \\
        | gzip > ${sample_id}_filtered.fastq.gz

    AFTER=\$(gunzip -c ${sample_id}_filtered.fastq.gz | grep -c '^@' || true)
    REMOVED=\$(( BEFORE - AFTER ))

    echo "Filtlong summary for ${sample_id}:"                          | tee -a ${sample_id}_filtlong.log
    echo "  Reads before : \${BEFORE}"                                 | tee -a ${sample_id}_filtlong.log
    echo "  Reads after  : \${AFTER} (min_length=${min_len} ${target_arg})" | tee -a ${sample_id}_filtlong.log
    echo "  Reads removed: \${REMOVED}"                                | tee -a ${sample_id}_filtlong.log

    [ -s ${sample_id}_filtered.fastq.gz ] || \\
        { echo "ERROR: No reads remain after filtering for ${sample_id}."; exit 1; }
    """

    stub:
    """
    gzip -c /dev/null > ${sample_id}_filtered.fastq.gz
    echo "stub log"   > ${sample_id}_filtlong.log
    """
}
