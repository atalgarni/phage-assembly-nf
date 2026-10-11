#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  ONT Phage Assembly Pipeline
//  Porechop_ABI → Filtlong → Flye → Medaka → Dnaapler → Stats / PhaBOX2 / Pharokka
//  Assembly, polishing, reorientation and annotation of phage genomes from ONT reads
//  Tuned for small genomes (5 kb – 200 kb) at high sequencing coverage
// ─────────────────────────────────────────────────────────────────────────────

include { PORECHOP_ABI    } from './modules/porechop_abi'
include { FILTLONG        } from './modules/filtlong'
include { FLYE            } from './modules/flye'
include { MEDAKA          } from './modules/medaka'
include { DNAAPLER        } from './modules/dnaapler'
include { PHABOX2         } from './modules/phabox2'
include { PHAROKKA        } from './modules/pharokka'
include { ASSEMBLY_STATS  } from './modules/assembly_stats'

// Optional host branch (--assemble_host)
include { HOST_READS        } from './modules/host_reads'
include { HOST_POOL         } from './modules/host_pool'
include { HOST_FLYE         } from './modules/host_flye'
include { HOST_MEDAKA       } from './modules/host_medaka'
include { KAPTIVE           } from './modules/kaptive'
include { MLST              } from './modules/mlst'
include { HOST_SUMMARY_ROW; HOST_CONCORDANCE } from './modules/host_summary'

// ── Parameter validation ──────────────────────────────────────────────────────

def validate_params() {
    if (!params.input) {
        error "ERROR: --input is required.\nUsage: nextflow run main.nf --input reads.fastq.gz --outdir results"
    }
    if (!params.phabox2_db)  error "ERROR: --phabox2_db is required"
    if (!params.pharokka_db) error "ERROR: --pharokka_db is required"

    def valid_read_types = ['nano-raw', 'nano-hq', 'nano-corr',
                            'pacbio-raw', 'pacbio-corr', 'pacbio-hifi']
    if (!valid_read_types.contains(params.read_type)) {
        error "ERROR: Invalid --read_type '${params.read_type}'.\nValid options: ${valid_read_types.join(', ')}"
    }
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {
    validate_params()

    // Print run summary to stdout
    log.info """
    ╔══════════════════════════════════════════════════╗
    ║  ONT Phage Assembly Pipeline                         ║
    ╚══════════════════════════════════════════════════╝
    Input reads   : ${params.input}
    Output dir    : ${params.outdir}
    Read type     : ${params.read_type}
    Genome size   : ${params.genome_size ?: 'auto-estimate (set --genome_size for phage!)'}
    ABI mode      : ${params.abi}
    No split      : ${params.no_split}   (recommended: true for phage)
    Medaka model  : ${params.medaka_model}
    PhaBOX2 DB    : ${params.phabox2_db}
    Pharokka DB   : ${params.pharokka_db}
    Host branch   : ${params.assemble_host ? "on (genome size ${params.host_genome_size}${params.pool_host ? ', + pooled' : ''})" : 'off (--assemble_host to enable)'}
    Threads       : ${params.threads}
    Profile       : ${workflow.profile}
    ══════════════════════════════════════════════════
    """.stripIndent()

    // Build input channel — supports single file or glob pattern
    // e.g. --input "reads/*.fastq.gz" will concatenate/treat as multiple samples
    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .map { file ->
            def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
            tuple(sample_id, file)
        }

    // ── Step 1: Adapter trimming ──────────────────────────────────────────────
    PORECHOP_ABI(ch_reads)

    // ── Step 2: Length filtering (removes reads < min_length for Flye) ──────────
    FILTLONG(PORECHOP_ABI.out.trimmed_reads)

    // ── Step 3: De novo assembly ──────────────────────────────────────────────
    FLYE(FILTLONG.out.filtered_reads)

    // ── Step 4: Polish with Medaka ────────────────────────────────────────────
    // Join FILTLONG filtered reads + FLYE assembly by sample_id
    // These are the exact same reads used for assembly
    ch_medaka = FILTLONG.out.filtered_reads
        .join(FLYE.out.assembly, by: 0)

    MEDAKA(ch_medaka)

    // ── Step 5: Reorient to start at terminase large subunit (terL) ─────────
    // Standardises orientation across all isolates before annotation.
    // Prevents reverse-assembly artifacts in pangenome and genome maps.
    DNAAPLER(MEDAKA.out.polished)

    // ── Steps 6-8 all take the REORIENTED FASTA and run in parallel ──────────
    // ── Step 6: Assembly statistics on reoriented assembly ────────────────────
    ASSEMBLY_STATS(DNAAPLER.out.reoriented)

    // ── Step 7: Phage characterisation (PhaBOX2) ──────────────────────────────
    PHABOX2(DNAAPLER.out.reoriented)

    // ── Step 8: Genome annotation (Pharokka) ──────────────────────────────────
    PHAROKKA(DNAAPLER.out.reoriented)

    // ── Steps 9–12 (optional): host genome from the same lysate reads ───────
    // Reads that do not map to the phage are host DNA carried over in the
    // lysate. They are assembled and typed (Kaptive K/O, MLST) per sample;
    // a concordance report shows which samples share a propagation host.
    if (params.assemble_host) {

        HOST_READS(PORECHOP_ABI.out.trimmed_reads.join(DNAAPLER.out.reoriented, by: 0))

        ch_host = HOST_READS.out.host_reads
            .filter { sample_id, reads, n ->
                def keep = (n as Integer) >= params.host_min_reads
                if (!keep) log.warn "Host branch: ${sample_id} has only ${n} host reads (< ${params.host_min_reads}) — skipped"
                keep
            }
            .map { sample_id, reads, n -> tuple(sample_id, reads) }

        if (params.pool_host) {
            HOST_POOL(HOST_READS.out.host_reads.map { it[1] }.collect())
            ch_host = ch_host.mix(HOST_POOL.out.host_reads)
        }

        HOST_FLYE(ch_host)
        HOST_MEDAKA(ch_host.join(HOST_FLYE.out.assembly, by: 0))

        ch_host_asm = HOST_MEDAKA.out.polished
        KAPTIVE(ch_host_asm)
        MLST(ch_host_asm)

        HOST_SUMMARY_ROW(
            ch_host_asm
                .join(KAPTIVE.out.results, by: 0)
                .join(MLST.out.results,    by: 0)
        )
        HOST_CONCORDANCE(HOST_SUMMARY_ROW.out.row.collect())
    }

    // ── Collect & report outputs ──────────────────────────────────────────────
    ASSEMBLY_STATS.out.stats
        .collectFile(name: 'assembly_summary.tsv',
                     storeDir: "${params.outdir}/summary",
                     keepHeader: true,
                     sort: { it[0] })
}

workflow.onComplete {
    log.info """
    ══════════════════════════════════════════════════
    Pipeline complete!
    Status    : ${workflow.success ? 'SUCCESS ✔' : 'FAILED ✘'}
    Duration  : ${workflow.duration}
    Results   : ${params.outdir}
    Report    : ${params.outdir}/pipeline_report.html
    ══════════════════════════════════════════════════
    """.stripIndent()
}
