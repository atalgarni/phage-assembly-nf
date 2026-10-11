// =============================================================================
//  Modules: HOST_SUMMARY_ROW + HOST_CONCORDANCE
//  One summary row per sample, then a concordance report across samples:
//  samples with identical KL / OL / ST profiles share the same host strain.
//  Plain shell + POSIX awk; no extra software needed.
// =============================================================================

process HOST_SUMMARY_ROW {
    tag "${sample_id}"
    label 'process_low'

    input:
    tuple val(sample_id), path(assembly), path(kaptive_k), path(kaptive_o), path(mlst)

    output:
    path "${sample_id}_host_row.tsv", emit: row

    script:
    """
    # Assembly metrics
    read CONTIGS TOTAL N50 <<< \$(awk '/^>/{if(s!=""){print length(s)}; s=""; next} {s=s\$0}
        END{if(s!=""){print length(s)}}' ${assembly} | sort -rn | \\
        awk '{l[NR]=\$1; t+=\$1} END{r=0; n50=0; for(i=1;i<=NR;i++){r+=l[i]; if(!n50 && r>=t/2) n50=l[i]}
             printf "%d %d %d", NR, t, n50}')

    # Kaptive: col2 = locus, col4 = confidence, col6 = identity
    KL=\$(awk -F'\\t' 'NR==2{print \$2}' ${kaptive_k});  KL=\${KL:--}
    KLC=\$(awk -F'\\t' 'NR==2{print \$4}' ${kaptive_k}); KLC=\${KLC:--}
    KLI=\$(awk -F'\\t' 'NR==2{print \$6}' ${kaptive_k}); KLI=\${KLI:--}
    OL=\$(awk -F'\\t' 'NR==2{print \$2}' ${kaptive_o});  OL=\${OL:--}
    OLC=\$(awk -F'\\t' 'NR==2{print \$4}' ${kaptive_o}); OLC=\${OLC:--}

    # mlst: col3 = ST, cols 4+ = alleles
    ST=\$(awk -F'\\t' 'END{print \$3}' ${mlst});        ST=\${ST:--}
    ALLELES=\$(awk -F'\\t' 'END{for(i=4;i<=NF;i++) printf "%s%s", \$i, (i<NF?",":"")}' ${mlst})

    echo -e "${sample_id}\\t\${CONTIGS}\\t\${TOTAL}\\t\${N50}\\t\${KL}\\t\${KLC}\\t\${KLI}\\t\${OL}\\t\${OLC}\\t\${ST}\\t\${ALLELES:--}" \\
        > ${sample_id}_host_row.tsv
    """

    stub:
    """
    echo -e "${sample_id}\\t12\\t5400000\\t900000\\tKL46\\tTypeable\\t99.5%\\tOL3\\tTypeable\\t3691\\t-" > ${sample_id}_host_row.tsv
    """
}

process HOST_CONCORDANCE {
    label 'process_low'

    publishDir(path: "${params.outdir}/09_host", mode: 'copy')

    input:
    path(rows)

    output:
    path "host_typing_summary.tsv", emit: summary
    path "host_concordance.txt",    emit: concordance

    script:
    """
    echo -e "sample\\tcontigs\\ttotal_bp\\tn50\\tKL_locus\\tKL_confidence\\tKL_identity\\tOL_locus\\tOL_confidence\\tMLST_ST\\tMLST_alleles" \\
        > host_typing_summary.tsv
    cat ${rows} | sort >> host_typing_summary.tsv

    # Group samples by capsule profile (KL + OL); STs listed per group.
    # Samples with an uncalled KL ("-") are reported as unresolved.
    awk -F'\\t' 'NR>1 && \$1!="pooled" {
            if (\$5 == "-" || \$5 == "") { unres = (unres == "" ? \$1 : unres ", " \$1); next }
            key = "KL=" \$5 "  OL=" \$8
            if (!(key in g)) { order[++k] = key }
            g[key] = (g[key] == "" ? \$1 : g[key] ", " \$1); n++
            if (\$10 != "-" && \$10 != "" && index(" " st[key] " ", " ST" \$10 " ") == 0)
                st[key] = (st[key] == "" ? "ST" \$10 : st[key] " ST" \$10)
        }
        END {
            print "Host concordance across samples"
            print "================================"
            for (i=1; i<=k; i++) {
                print "Profile " i ": " order[i] "  MLST: " (st[order[i]] == "" ? "not called" : st[order[i]])
                print "  Samples: " g[order[i]]
            }
            if (unres != "") print "Unresolved (capsule locus not typed): " unres
            print ""
            if (k == 1 && n > 1)
                print "All typed samples share one host profile -> consistent with a single propagation host.\\nPooling host reads (--pool_host) can improve host assembly completeness."
            else if (k > 1)
                print k " distinct host profiles -> samples were propagated on different hosts; do not pool across profiles."
            print "Uncalled values (-) usually reflect low host depth and are not evidence of a different host."
        }' host_typing_summary.tsv > host_concordance.txt

    cat host_concordance.txt
    """

    stub:
    """
    cat ${rows} > host_typing_summary.tsv
    echo stub > host_concordance.txt
    """
}
