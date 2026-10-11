// =============================================================================
//  Module: HOST_FLYE
//  De novo assembly of host reads (bacterial genome size, e.g. 5.5m).
//  A failed or fragmented host assembly never stops the phage pipeline.
// =============================================================================

process HOST_FLYE {
    tag "${sample_id}"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/09_host/${sample_id}/flye",
        mode: 'copy',
        saveAs: { filename -> filename.endsWith('.log') ? "logs/${filename}" : filename }
    )

    input:
    tuple val(sample_id), path(host_reads)

    output:
    tuple val(sample_id), path("${sample_id}_host_assembly.fasta"), emit: assembly
    path "${sample_id}_host_assembly_info.txt",                      emit: info
    path "${sample_id}_host_flye.log",                               emit: log

    script:
    """
    flye \\
        --${params.read_type} ${host_reads} \\
        --genome-size ${params.host_genome_size} \\
        --out-dir     flye_host \\
        --threads     ${task.cpus} \\
        2>&1 | tee ${sample_id}_host_flye.log

    cp flye_host/assembly.fasta     ${sample_id}_host_assembly.fasta
    cp flye_host/assembly_info.txt  ${sample_id}_host_assembly_info.txt

    [ -s ${sample_id}_host_assembly.fasta ] || \\
        { echo "ERROR: Flye produced no host assembly for ${sample_id}"; exit 1; }
    """

    stub:
    """
    printf '>contig_1\\nACGTACGT\\n' > ${sample_id}_host_assembly.fasta
    touch ${sample_id}_host_assembly_info.txt ${sample_id}_host_flye.log
    """
}
