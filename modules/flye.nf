// =============================================================================
//  Module: FLYE
//  De novo assembly of trimmed ONT reads
// =============================================================================

process FLYE {
    tag "${sample_id}"
    label 'process_high'

    publishDir(
        path: "${params.outdir}/03_flye_assembly/${sample_id}",
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log'))          "logs/${filename}"
            else if (filename == 'assembly.fasta')  filename
            else if (filename == 'assembly_graph.gfa') filename
            else if (filename == 'assembly_info.txt')  filename
            else                                    filename
        }
    )

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("assembly.fasta"),        emit: assembly
    path "assembly_graph.gfa",                           emit: graph
    path "assembly_info.txt",                            emit: info
    path "${sample_id}_flye.log",                        emit: log

    script:
    def genome_size_arg = params.genome_size  ? "--genome-size ${params.genome_size}" : ''
    def meta_flag       = params.meta         ? '--meta'                              : ''
    def scaffold_flag   = params.scaffold     ? '--scaffold'                          : ''
    def min_overlap_arg = params.min_overlap  ? "--min-overlap ${params.min_overlap}" : ''

    """
    flye \\
        --${params.read_type} ${trimmed_reads} \\
        --out-dir    flye_out \\
        --threads    ${task.cpus} \\
        --iterations ${params.iterations} \\
        ${genome_size_arg} \\
        ${meta_flag} \\
        ${scaffold_flag} \\
        ${min_overlap_arg} \\
        2>&1 | tee ${sample_id}_flye.log

    # Move key outputs to process working directory for publishDir
    cp flye_out/assembly.fasta        assembly.fasta
    cp flye_out/assembly_graph.gfa    assembly_graph.gfa
    cp flye_out/assembly_info.txt     assembly_info.txt

    # Fail loudly if assembly is empty
    [ -s assembly.fasta ] || \\
        { echo "ERROR: Flye produced no assembly for ${sample_id}"; exit 1; }
    """

    stub:
    """
    echo "Stub: Flye for ${sample_id}"
    echo ">contig_1" > assembly.fasta
    echo "ATCGATCGATCG" >> assembly.fasta
    touch assembly_graph.gfa
    echo -e "seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path" > assembly_info.txt
    echo "stub log" > ${sample_id}_flye.log
    """
}
