// =============================================================================
//  Module: ASSEMBLY_STATS
//  Compute basic assembly metrics: contigs, total bp, N50, N90, GC%
//  Built-in: POSIX awk only, no external tool (identical output on every system).
//  POSIX awk only — no asort(), no gawk required (macOS BSD awk compatible)
// =============================================================================

process ASSEMBLY_STATS {
    tag "${sample_id}"
    label 'process_low'

    publishDir(
        path: "${params.outdir}/06_assembly_stats",
        mode: 'copy'
    )

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}_assembly_stats.tsv", emit: stats

    script:
    """
    # ── Pure POSIX awk + sort (no asort / no gawk needed) ─────────────────
    # Pass 1: extract per-contig length and GC, write lengths to temp file
    awk '
    /^>/ {
        if (seq != "") {
            len = length(seq)
            total += len
            gc_tmp = seq
            gsub(/[^GCgc]/, "", gc_tmp)
            gc += length(gc_tmp)
            print len > "contig_lengths.txt"
            n++
        }
        seq = ""; next
    }
    { seq = seq \$0 }
    END {
        if (seq != "") {
            len = length(seq)
            total += len
            gc_tmp = seq
            gsub(/[^GCgc]/, "", gc_tmp)
            gc += length(gc_tmp)
            print len > "contig_lengths.txt"
            n++
        }
        # Write totals for pass 2
        print total > "total_bp.txt"
        print gc    > "total_gc.txt"
        print n     > "total_n.txt"
    }
    ' ${assembly}

    TOTAL_BP=\$(cat total_bp.txt)
    TOTAL_GC=\$(cat total_gc.txt)
    TOTAL_N=\$(cat total_n.txt)

    # Pass 2: sort lengths descending, compute N50 / N90 / largest
    echo -e "sample\\tnum_contigs\\ttotal_bp\\tlargest_contig\\tn50\\tn90\\tgc_pct" \\
        > ${sample_id}_assembly_stats.tsv

    sort -rn contig_lengths.txt | awk \\
        -v sample="${sample_id}" \\
        -v total="\${TOTAL_BP}" \\
        -v gc="\${TOTAL_GC}" \\
        -v nctg="\${TOTAL_N}" \\
    '
    BEGIN { run=0; n50=0; n90=0; largest=0; first=1 }
    {
        if (first) { largest=\$1; first=0 }
        run += \$1
        if (!n50 && run >= total*0.5) n50=\$1
        if (!n90 && run >= total*0.9) n90=\$1
    }
    END {
        gc_pct = (total>0) ? gc/total*100 : 0
        printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%.2f\\n",
               sample, nctg, total, largest, n50, n90, gc_pct
    }' >> ${sample_id}_assembly_stats.tsv

    rm -f contig_lengths.txt total_bp.txt total_gc.txt total_n.txt

    echo "────────────────────────────────────────"
    echo "Assembly stats for ${sample_id}:"
    cat ${sample_id}_assembly_stats.tsv
    echo "────────────────────────────────────────"
    """

    stub:
    """
    echo -e "sample\\tnum_contigs\\ttotal_bp\\tlargest_contig\\tn50\\tn90\\tgc_pct" \\
        > ${sample_id}_assembly_stats.tsv
    echo -e "${sample_id}\\t1\\t110000\\t110000\\t110000\\t110000\\t45.00" \\
        >> ${sample_id}_assembly_stats.tsv
    """
}
