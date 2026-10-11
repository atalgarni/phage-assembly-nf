# ONT Phage Assembly Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-brightgreen)
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

A Nextflow DSL2 pipeline for de novo assembly, polishing, reorientation and annotation of bacteriophage genomes from Oxford Nanopore reads, tuned for single phage isolates (5–200 kb) sequenced at high coverage on R10.4.1 flow cells.

**From one lysate run, two genomes.** Host DNA carried over in the phage lysate is usually discarded as contamination. With `--assemble_host`, reads that do not map to the phage are assembled into a host genome and typed (capsule K/O locus with Kaptive, MLST), and a concordance report shows which phages were propagated on the same host strain — without sequencing the host separately.

```
Raw ONT reads (.fastq.gz)
        │
        ▼
 PORECHOP_ABI    ab initio adapter discovery & trimming (--no_split by default)
        ▼
 FILTLONG        length filter (≥1 kb) + optional longest-first subsampling
        ▼
 FLYE            de novo assembly (no --meta; circularity reported)
        ▼
 MEDAKA          consensus polishing (R10.4.1 model)
        ▼
 DNAAPLER        reorient genome to start at terL
        │
        ├──► ASSEMBLY_STATS   contigs, length, N50/N90, GC%
        ├──► PHABOX2          taxonomy, lifestyle, host prediction
        ├──► PHAROKKA         functional annotation, GenBank/GFF
        │
        └──► (--assemble_host)
             HOST_READS       trimmed reads not mapping to the phage (minimap2 + samtools)
             HOST_FLYE        host assembly (e.g. 5.5 Mb)  → HOST_MEDAKA polishing
             KAPTIVE + MLST   capsule K/O locus and sequence type
             HOST_CONCORDANCE which samples share a host profile
             (--pool_host)    extra assembly from pooled host reads of all samples
```

## Design choices for phage data

| Choice | Default | Why |
|---|---|---|
| `--no_split` | on | Chimera splitting on a small genome risks fragmenting legitimate long reads |
| `--meta` | off | A phage isolate is a single genome, not a community |
| `--genome_size` | *set it* | Flye's auto-estimate is less reliable at the extreme depths typical of phage runs |
| `--target_bases` | off | Subsample very deep or host-contaminated samples (e.g. `11000000` ≈ 100× of a 110 kb phage) |
| Dnaapler | on | Consistent start coordinate (terL) across isolates for comparison and annotation |

## Requirements

- Nextflow ≥ 24.04
- Conda/mamba (recommended), Docker (containers built via Seqera Wave), or all tools in `PATH`
- Dnaapler can share the Pharokka environment
- Databases (required): [PhaBOX2](https://github.com/KennthShang/PhaBOX) and [Pharokka](https://github.com/gbouras13/pharokka)

## Quick start

```bash
# Check the pipeline logic without running any tool
nextflow run main.nf -profile test -stub

# Full run, conda environments built automatically from bioconda
nextflow run main.nf \
    --input "data/*.fastq.gz" \
    --genome_size 50k \
    --phabox2_db /path/to/phabox_db \
    --pharokka_db /path/to/pharokka_db \
    --threads 16 \
    -profile conda
```

Each FASTQ file is treated as one sample; samples run in parallel.

### Using existing conda environments (`run.sh`)

If the tools are already installed in named conda environments, edit the CONFIG section at the top of `run.sh` (input, genome size, environment names, database paths) and run:

```bash
bash run.sh
```

`run.sh` checks every tool and input before launching, then reports the circularity of each assembly.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | *(required)* | Reads file or glob, e.g. `"data/*.fastq.gz"` |
| `--outdir` | `results` | Output directory |
| `--genome_size` | `null` | Expected phage genome size, e.g. `50k`, `150k` |
| `--read_type` | `nano-raw` | Flye read type without leading dashes: `nano-raw`, `nano-hq` (R10.4.1 HAC/SUP, Q20+), `nano-corr` |
| `--min_length` | `1000` | Filtlong minimum read length (bp) |
| `--target_bases` | `null` | Filtlong subsampling target (bases) |
| `--iterations` | `1` | Flye polishing iterations |
| `--min_overlap` | `null` | Flye minimum overlap (auto if unset) |
| `--medaka_model` | `r1041_e82_400bps_hac_v4.3.0` | Must match the basecalling model |
| `--phabox2_db` | *(required)* | PhaBOX2 database directory |
| `--pharokka_db` | *(required)* | Pharokka database directory |
| `--abi` | `true` | Porechop_ABI ab initio adapter discovery |
| `--no_split` | `true` | Disable chimera splitting |
| `--phabox2_skip_phamer` | `Y` | Skip PhaMer virus identification (isolates are known phages) |
| `--threads` | `8` | CPUs per process |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `16.GB` / `48.h` | Upper limits for any single step; raise on larger machines |

### Host branch

| Parameter | Default | Description |
|---|---|---|
| `--assemble_host` | `false` | Assemble and type the host from reads not mapping to the phage |
| `--pool_host` | `false` | Also assemble the pooled host reads of all samples (only when they share one host) |
| `--host_genome_size` | `5.5m` | Expected host genome size |
| `--host_min_reads` | `1000` | Skip host assembly for samples with fewer host reads |
| `--kaptive_k_db` / `--kaptive_o_db` | `kpsc_k` / `kpsc_o` | Kaptive databases (`ab_k` / `ab_o` for *A. baumannii*) |
| `--mlst_scheme` | `klebsiella` | MLST scheme; `null` for automatic detection |

## Output

```
results/
├── 01_porechop_abi/          trimmed reads + logs
├── 02_filtlong/              filtered reads + logs
├── 03_flye_assembly/<sample>/
│   ├── assembly.fasta        draft assembly
│   ├── assembly_graph.gfa    view in Bandage
│   └── assembly_info.txt     per-contig stats incl. circularity (circ. = Y)
├── 04_medaka/<sample>/       polished genome
├── 05_dnaapler/              final reoriented genome  ← use this
├── 06_assembly_stats/        per-sample statistics
├── 07_phabox2/<sample>/      taxonomy, lifestyle, host prediction
├── 08_pharokka/<sample>/     annotation (.gbk, .gff, .faa)
├── 09_host/                  (--assemble_host)
│   ├── <sample>/             host reads + depth, assembly, Kaptive, MLST
│   ├── pooled/               (--pool_host)
│   ├── host_typing_summary.tsv
│   └── host_concordance.txt
├── summary/assembly_summary.tsv
└── pipeline_report.html, pipeline_timeline.html, pipeline_dag.html, pipeline_trace.txt
```

A single contig close to the expected genome size with high, even coverage indicates a complete phage genome. Do not rely on circularity alone: many phages have linear genomes with terminal repeats or cohesive ends, so check `assembly_info.txt` together with comparison to the closest reference genome.

### Host branch: what to expect

Host depth depends on how much host DNA remains in the lysate (DNase treatment removes most of it), so host assemblies are often fragmented. Capsule typing still works on fragmented assemblies because Kaptive only needs the ~30 kb capsule locus; MLST and gene-level conclusions need more depth. Check `expected_depth` in `09_host/<sample>/<sample>_host_reads.stats.tsv`. Uncalled values (`-`) usually reflect low depth, not a different host. When all samples share one profile, `--pool_host` combines their host reads for a more complete host genome.

## Resume a run

```bash
nextflow run main.nf --input "data/*.fastq.gz" --genome_size 50k -profile conda -resume
```

## Tested on

Nextflow 25.10.0 on macOS (Intel) <!-- TO FILL: macOS version, e.g. 14.5 -->

| Tool | Version | | Tool | Version |
|---|---|---|---|---|
| Porechop_ABI | 0.5.0 | | Pharokka | 1.7.4 |
| Filtlong | 0.3.1 | | minimap2 | 2.28 |
| Flye | 2.9.6 | | SAMtools | 1.20 |
| Medaka | 2.1.1 | | Kaptive | 3.1.0 |
| Dnaapler | 1.4.0 | | mlst | 2.33.1 |
| PhaBOX2 | 2.1.13 | | Assembly statistics | built-in (awk) |

## Citation

If you use this pipeline, please cite it using the DOI above (see `CITATION.cff`), together with the tools it wraps: Porechop_ABI, Filtlong, Flye, Medaka, Dnaapler, PhaBOX2, Pharokka and, for the host branch, minimap2, SAMtools, Kaptive and mlst.

## Author

Abdulmalek T. Algarni — King Abdulaziz City for Science and Technology (KACST)

## License

MIT — see `LICENSE`.
