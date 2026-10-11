// =============================================================================
//  Module: HOST_READS
//  Extract host reads = trimmed reads that do NOT map to the final phage genome
//  (minimap2 map-ont → samtools -f 4). Reports read count, bases and the
//  expected host coverage for the configured host genome size.
// =============================================================================

process HOST_READS {
    tag "${sample_id}"
    label 'process_medium'

    publishDir(
        path: "${params.outdir}/09_host/${sample_id}",
        mode: 'copy'
    )

    input:
    tuple val(sample_id), path(reads), path(phage)

    output:
    tuple val(sample_id), path("${sample_id}_host_reads.fastq.gz"), env(HOST_N), emit: host_reads
    path "${sample_id}_host_reads.stats.tsv",                                   emit: stats

    script:
    def gsize = params.host_genome_size
    """
    minimap2 -ax map-ont -t ${task.cpus} ${phage} ${reads} 2> /dev/null \\
        | samtools view -f 4 -b \\
        | samtools fastq \\
        | gzip > ${sample_id}_host_reads.fastq.gz

    # Read count and total bases (POSIX awk)
    gunzip -c ${sample_id}_host_reads.fastq.gz \\
        | awk 'NR%4==2 {n++; b+=length(\$0)} END {printf "%d\\t%d\\n", n+0, b+0}' > counts.txt
    HOST_N=\$(cut -f1 counts.txt)
    HOST_BP=\$(cut -f2 counts.txt)

    # Expected depth = bases / genome size ("5.5m" → 5500000)
    GSIZE=\$(echo "${gsize}" | awk '{v=\$0; m=1; u=tolower(substr(v,length(v)))
             if(u=="k"){m=1e3; v=substr(v,1,length(v)-1)}
             else if(u=="m"){m=1e6; v=substr(v,1,length(v)-1)}
             else if(u=="g"){m=1e9; v=substr(v,1,length(v)-1)}
             printf "%d", v*m}')
    DEPTH=\$(awk -v b=\${HOST_BP} -v g=\${GSIZE} 'BEGIN{printf "%.1f", (g>0)? b/g : 0}')

    echo -e "sample\\thost_reads\\thost_bases\\texpected_depth" > ${sample_id}_host_reads.stats.tsv
    echo -e "${sample_id}\\t\${HOST_N}\\t\${HOST_BP}\\t\${DEPTH}" >> ${sample_id}_host_reads.stats.tsv
    echo "Host reads for ${sample_id}: \${HOST_N} reads, \${HOST_BP} bp (~\${DEPTH}x of ${gsize})"
    """

    stub:
    """
    printf '@r1\\nACGT\\n+\\nIIII\\n' | gzip > ${sample_id}_host_reads.fastq.gz
    HOST_N=5000
    echo -e "sample\\thost_reads\\thost_bases\\texpected_depth\\n${sample_id}\\t5000\\t50000000\\t9.1" > ${sample_id}_host_reads.stats.tsv
    """
}
