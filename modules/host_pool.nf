// =============================================================================
//  Module: HOST_POOL
//  Concatenate host reads from all samples into one "pooled" sample.
//  Use only when the samples share the same propagation host
//  (check 09_host/host_concordance.txt from a per-sample run first).
// =============================================================================

process HOST_POOL {
    tag "pooled"
    label 'process_low'

    input:
    path(host_reads)

    output:
    tuple val("pooled"), path("pooled_host_reads.fastq.gz"), emit: host_reads

    script:
    """
    cat ${host_reads} > pooled_host_reads.fastq.gz
    """

    stub:
    """
    cat ${host_reads} > pooled_host_reads.fastq.gz
    """
}
