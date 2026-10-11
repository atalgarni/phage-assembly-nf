# phage-assembly-nf

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23300714.svg)](https://doi.org/10.5281/zenodo.23300714)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-brightgreen)
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

A Nextflow DSL2 workflow for long-read-only reconstruction and characterisation of bacteriophage genomes from Oxford Nanopore Technologies (ONT) sequencing of phage lysates, with an optional branch that recovers and types the **propagation host** from the same run.

For each isolate the workflow performs adapter removal, length-based read selection, de novo assembly, neural-network consensus polishing, reorientation to the large terminase subunit (*terL*), assembly statistics, taxonomic/lifestyle/host prediction and functional annotation. Reads that do not map to the final phage genome — residual host DNA that survives DNase treatment of the lysate — can optionally be assembled into a draft host genome and typed by capsule (K) and O-antigen (O) locus and multilocus sequence type (MLST). A cross-sample concordance report identifies isolates propagated on the same host strain.

---

## Contents

1. [Workflow overview](#1-workflow-overview)
2. [Methods](#2-methods)
3. [Host branch](#3-host-branch)
4. [Installation](#4-installation)
5. [Input](#5-input)
6. [Usage](#6-usage)
7. [Parameters](#7-parameters)
8. [Output](#8-output)
9. [Interpreting the results](#9-interpreting-the-results)
10. [Validation](#10-validation)
11. [Computational requirements](#11-computational-requirements)
12. [Troubleshooting](#12-troubleshooting)
13. [Software versions](#13-software-versions)
14. [Citation](#14-citation)
15. [References](#15-references)

---

## 1. Workflow overview

```
                         ONT reads, one FASTQ per isolate
                                       │
                       ┌───────────────▼───────────────┐
                       │ PORECHOP_ABI   adapter removal │
                       └───────────────┬───────────────┘
                                       │ trimmed reads ─────────────────────┐
                       ┌───────────────▼───────────────┐                    │
                       │ FILTLONG  length filter /      │                    │
                       │           depth subsampling    │                    │
                       └───────────────┬───────────────┘                    │
                                       │ filtered reads ──────┐              │
                       ┌───────────────▼───────────────┐      │              │
                       │ FLYE   de novo assembly        │      │              │
                       └───────────────┬───────────────┘      │              │
                       ┌───────────────▼───────────────┐      │              │
                       │ MEDAKA consensus polishing     │◄─────┘              │
                       └───────────────┬───────────────┘                    │
                       ┌───────────────▼───────────────┐                    │
                       │ DNAAPLER reorient to terL      │                    │
                       └───────────────┬───────────────┘                    │
                                       │ final phage genome                  │
          ┌────────────────────────────┼────────────────────────┐           │
          ▼                            ▼                        ▼           ▼
   ASSEMBLY_STATS               PHABOX2                     PHAROKKA     HOST_READS (--assemble_host)
   length, N50/N90, GC     taxonomy (PhaGCN)            CDS calling      minimap2 map-ont → unmapped reads
                           lifestyle (PhaTYP)           PHROG/VFDB/CARD          │
                           host (CHERRY)                tRNA/CRISPR/tmRNA  HOST_FLYE → HOST_MEDAKA
                                                        INPHARED Mash            │
                                                                          KAPTIVE (K/O)  +  MLST
                                                                                 │
                                                                    HOST_SUMMARY + HOST_CONCORDANCE
```

Samples are processed independently and in parallel; each FASTQ file is one sample.

---

## 2. Methods

### 2.1 Adapter removal — Porechop_ABI 0.5.0 [2]

```
porechop_abi -abi -i <reads> -o <sample>_trimmed.fastq.gz --format fastq.gz \
             --no_split --min_split_read_size 1000 -t <cpus>
```

* `-abi` infers adapter sequences *ab initio* from the reads rather than relying only on Porechop's built-in adapter set, which is useful for newer ligation and rapid kits.
* `--no_split` (default) disables splitting of reads at internal adapters. Phage genomes are small and sequenced at very high depth, so long reads spanning most of the genome are the most informative for assembly; splitting would fragment legitimate long reads in exchange for removing a small number of chimeras. Set `--no_split false` to restore chimera splitting.

### 2.2 Read selection — Filtlong 0.3.1 [3]

```
filtlong --min_length 1000 [--target_bases N] <trimmed reads>
```

* Reads shorter than `--min_length` (default 1 kb) are discarded. This matches Flye's minimum overlap floor; shorter reads cannot contribute to the overlap graph.
* `--target_bases` (optional) retains the highest-scoring reads (Filtlong ranks by length and mean quality) up to the requested number of bases. Phage lysates frequently yield 1,000–10,000× coverage; downsampling to ~100× (e.g. `--target_bases 11000000` for a ~110 kb genome) reduces run time and avoids coverage-driven graph artefacts.
* Read counts before and after filtering are written to `02_filtlong/<sample>_filtlong.log`.

### 2.3 De novo assembly — Flye 2.9.6 [4]

```
flye --nano-raw <filtered reads> --genome-size <size> --iterations 1 \
     --out-dir flye_out --threads <cpus>
```

* Run in single-genome mode (no `--meta`): an isolate lysate contains one phage genome, and metagenomic mode would relax the uniform-coverage assumption that helps resolve a single high-depth genome.
* `--genome-size` should be set to the expected phage genome size. Flye uses it to estimate coverage and choose the initial disjointig parameters; its automatic estimate is less reliable at the extreme depths typical of phage runs.
* `--read_type` sets the Flye input mode. The default is `nano-raw`; `nano-hq` is recommended by the Flye authors for R10.4.1 data basecalled in HAC/SUP mode (Q20+).
* Flye is retried up to twice with increasing memory (8 GB × attempt).

### 2.4 Consensus polishing — Medaka 2.1.1 [5]

```
medaka_consensus -i <filtered reads> -d <flye assembly> -m r1041_e82_400bps_hac_v4.3.0 \
                 -o medaka_out -t <cpus>
```

* Filtered reads are realigned to the Flye draft and a recurrent neural network computes a consensus, correcting the systematic indel errors of nanopore reads (predominantly in homopolymers).
* **The model must match the flow cell, kit and basecaller model used to generate the reads.** The default corresponds to Dorado `dna_r10.4.1_e8.2_400bps_hac@v4.3.0`. List available models with `medaka tools list_models`.

### 2.5 Reorientation — Dnaapler 1.4.0 [6]

```
dnaapler phage -i <polished assembly> -o dnaapler_out -p <sample> -t <cpus>
```

* Rotates the assembly so that it begins at the start codon of the large terminase subunit (*terL*), identified by MMseqs2 search against a *terL* database, and places that gene on the forward strand.
* This gives every isolate the same coordinate origin, which is required for meaningful gene-order (synteny) comparison, pangenome analysis and consistent locus numbering in annotations.
* If Dnaapler fails or finds no *terL* hit, the polished assembly is passed on unchanged and a warning is written to the log; the pipeline does not stop.

### 2.6 Assembly statistics — built-in

Number of contigs, total length, largest contig, N50, N90 and GC content are computed with POSIX `awk` (no external dependency, identical output on Linux and macOS). Per-sample tables are concatenated into `summary/assembly_summary.tsv`.

### 2.7 Taxonomy, lifestyle and host prediction — PhaBOX2 2.1.13 [7–10]

```
phabox2 --task end_to_end --contigs <genome> --dbdir <db> --len 1000 --skip Y
```

| Module | Prediction | Basis |
|---|---|---|
| PhaGCN [8] | Taxonomic lineage to genus (and nearest species) | Graph convolutional network over protein-sharing networks with ICTV reference genomes |
| PhaTYP [9] | Lifestyle (virulent / temperate) | Transformer model on protein-cluster "sentences" |
| CHERRY [10] | Host genus | Link prediction on a phage–host knowledge graph, with alignment-based (CRISPR, AAI) evidence where available |

`--skip Y` bypasses PhaMer (virus identification), because input contigs are assembled from isolate lysates and are already known to be phage. pandas is pinned to 2.2 in the PhaBOX2 environment: PhaBOX2 2.1.13 raises a `TypeError` in PhaGCN under pandas 3.

### 2.8 Functional annotation — Pharokka 1.10.1 [11]

```
pharokka run -i <genome> -o pharokka_out -d <db> -g phanotate -p <sample> -t <cpus> -f
```

* **Gene calling:** PHANOTATE [12], which models phage gene density and permits overlapping genes, and generally calls more short ORFs than Prodigal on phage genomes.
* **Functional assignment:** MMseqs2 [13] profile search and PyHMMER [14] HMM search against PHROGs [15] (functional categories: head and packaging, tail, lysis, DNA/RNA metabolism, integration/excision, transcription regulation, moron/auxiliary metabolic genes, connector, other, unknown).
* **Additional databases:** VFDB [16] (virulence factors) and CARD [17] (antimicrobial resistance genes) — relevant to therapeutic phage safety assessment.
* **Non-coding features:** tRNAscan-SE 2 [18] (tRNA), Aragorn [19] (tmRNA), MinCED [20] (CRISPR arrays).
* **Closest known phage:** Mash [21] search against INPHARED [22] (all complete phage genomes in GenBank).
* The wrapper detects the CLI generation: `pharokka run` (≥ 1.10) or `pharokka.py` (older).

---

## 3. Host branch

Enabled with `--assemble_host`. Designed for lysates in which host DNA was not completely removed.

### 3.1 Host read extraction — minimap2 [23], SAMtools [24]

```
minimap2 -ax map-ont <final phage genome> <Porechop-trimmed reads> \
  | samtools view -f 4 -b | samtools fastq | gzip > <sample>_host_reads.fastq.gz
```

* Reads are taken **after Porechop_ABI and before Filtlong**, so short host fragments that Filtlong would discard are retained.
* `-f 4` keeps reads that do not align to the phage. These are predominantly host chromosomal and plasmid DNA; contaminants and unassembled phage reads, if present, are also included.
* Expected host depth is computed as `host bases ÷ --host_genome_size` and written to `<sample>_host_reads.stats.tsv`.

### 3.2 Depth criteria

| Condition | Behaviour |
|---|---|
| host reads < `--host_min_reads` (1,000) | Host assembly skipped, warning logged |
| expected depth < `--host_min_depth` (10×) | Host **is** assembled, warning logged, sample flagged `LOW` in the summary |
| expected depth ≥ `--host_min_depth` | Assembled and flagged `OK` |

Low-depth samples are not skipped because a fragmented host assembly still supports species-level identification (e.g. by BLAST of the largest contigs), whereas K/O-locus and MLST calls become unreliable when the relevant loci are incompletely assembled.

### 3.3 Assembly, polishing and typing

* **Flye** in single-genome mode with `--genome-size <host_genome_size>` (default 5.5 Mb, *K. pneumoniae*), then **Medaka** with the same model as the phage. If Medaka fails, the Flye draft is used.
* **Kaptive 3.1.0** [25] identifies the capsule polysaccharide (K) and lipopolysaccharide O-antigen (O) biosynthesis loci and reports the best-matching locus, match confidence (e.g. *Typeable*) and identity. Default databases are `kpsc_k` / `kpsc_o` (*Klebsiella pneumoniae* species complex); use `ab_k` / `ab_o` for *Acinetobacter baumannii*.
* **mlst 2.33.1** [26] (PubMLST schemes [27]) assigns the sequence type; default scheme `klebsiella`, or automatic scheme detection with `--mlst_scheme false`.
* Failure of any host step affects only that sample and never stops the phage workflow.

### 3.4 Concordance report

`HOST_CONCORDANCE` groups samples by capsule profile (KL + OL locus) and lists the STs observed in each group:

* **One profile across samples** → consistent with a single propagation host; `--pool_host` can then improve completeness.
* **Several profiles** → samples were propagated on different hosts; host reads should not be pooled across profiles.
* Samples without a K-locus call are reported as *unresolved*, and low-depth samples are listed separately. An uncalled value (`-`) usually reflects low depth and is not evidence of a different host.

### 3.5 Pooled host assembly

With `--pool_host`, host reads from **all** samples are concatenated and assembled as an additional sample named `pooled`, increasing depth when the isolates share a host. Use only after the per-sample concordance report supports a common host.

---

## 4. Installation

### 4.1 Nextflow

Nextflow [1] ≥ 24.04 and Java 17–26 are required. A dedicated environment avoids conflicts with an older Java in the base conda environment:

```bash
mamba create -n nextflow -c conda-forge -c bioconda nextflow=25.10.0 openjdk=17 -y
conda activate nextflow
nextflow -version
```

### 4.2 Pipeline

```bash
git clone https://github.com/atalgarni/phage-assembly-nf.git
cd phage-assembly-nf
nextflow run main.nf -profile test -stub --assemble_host   # logic check, no tools executed
```

### 4.3 Software

| Profile | How tools are provided | Status |
|---|---|---|
| `conda` | Nextflow builds one environment per process from the pinned Bioconda [28] versions (section 13). Requires `mamba`. | Tested (macOS Intel) |
| `docker` | Containers built on the fly from the same conda specifications via Seqera Wave (internet access required) | Not yet tested |
| `standard` | All tools must already be in `PATH` | — |

Existing environments can be used instead of building new ones by passing full environment paths (`--flye_conda /path/to/envs/flye`, etc.); `run.sh` does this automatically from environment names.

### 4.4 Databases

**PhaBOX2** — download the database matching PhaBOX2 2.1.x from the [PhaBOX repository](https://github.com/KennthShang/PhaBOX) and pass its directory with `--phabox2_db`.

**Pharokka** — the database must match the Pharokka release used by the pipeline (1.10.1). A database from another release fails Pharokka's version check (e.g. `PHROGs Database file VERSION_1_8_0 is missing`). Install into a new, empty folder:

```bash
mamba create -n pharokka110 -c conda-forge -c bioconda pharokka=1.10.1 -y
mamba run -n pharokka110 pharokka install -o pharokka_db_v110
```

Both database directories are checked for existence and content at launch, before any assembly is run.

**Kaptive and mlst** databases are bundled with their conda packages.

---

## 5. Input

* One FASTQ (`.fastq`, `.fq`, optionally gzipped) per isolate, passed as a path or a quoted glob: `--input "reads/*.fastq.gz"`.
* The sample ID is the file name without the FASTQ extension (`T_7.1_A.fastq.gz` → `T_7.1_A`).
* Reads should be basecalled and demultiplexed. Concatenate multiple FASTQ chunks of one barcode beforehand (`cat barcode01/*.fastq.gz > sample.fastq.gz`).
* Note the basecalling model: it determines `--medaka_model` and the appropriate `--read_type`.

---

## 6. Usage

**Single isolate, phage only**

```bash
nextflow run main.nf \
    --input reads/T_7.1_A.fastq.gz \
    --genome_size 111k \
    --phabox2_db /db/phabox_db_v2 \
    --pharokka_db /db/pharokka_db_v110 \
    --outdir results \
    -profile conda
```

**Several isolates with host recovery and typing**

```bash
nextflow run main.nf \
    --input "reads/*.fastq.gz" \
    --genome_size 111k \
    --target_bases 11000000 \
    --phabox2_db /db/phabox_db_v2 \
    --pharokka_db /db/pharokka_db_v110 \
    --assemble_host \
    --outdir results \
    -profile conda
```

**Pooled host assembly** (after the concordance report shows a single host profile): add `--pool_host`.

**Non-*Klebsiella* host** (e.g. *A. baumannii*, ~4 Mb):

```bash
    --assemble_host --host_genome_size 4m --kaptive_k_db ab_k --kaptive_o_db ab_o --mlst_scheme abaumannii_2
```

**Resume** an interrupted or modified run (completed steps are reused from `work/`):

```bash
nextflow run main.nf <same arguments> -resume
```

**Wrapper script:** edit the CONFIG block at the top of `run.sh` (input, genome size, environment names, database paths) and run `bash run.sh`. It checks every tool and input before launching and prints assembly topology at the end.

---

## 7. Parameters

### Input / output

| Parameter | Default | Description |
|---|---|---|
| `--input` | *(required)* | FASTQ file or quoted glob |
| `--outdir` | `results` | Output directory |

### Read processing

| Parameter | Default | Description |
|---|---|---|
| `--abi` | `true` | Porechop_ABI *ab initio* adapter inference |
| `--no_split` | `true` | Disable splitting of reads at internal adapters |
| `--discard_middle` | `false` | Discard reads with internal adapters (when splitting is enabled) |
| `--min_split_read_size` | `1000` | Minimum length of split read fragments |
| `--min_length` | `1000` | Filtlong minimum read length (bp) |
| `--target_bases` | `null` | Filtlong target bases (off by default) |

### Assembly and polishing

| Parameter | Default | Description |
|---|---|---|
| `--genome_size` | `null` | Expected phage genome size, e.g. `50k`, `111k` (strongly recommended) |
| `--read_type` | `nano-raw` | Flye input mode without leading dashes: `nano-raw`, `nano-hq`, `nano-corr` |
| `--iterations` | `1` | Flye polishing iterations |
| `--min_overlap` | `null` | Flye minimum overlap (automatic if unset) |
| `--meta` | `false` | Flye metagenome mode (not recommended for isolates) |
| `--scaffold` | `false` | Flye scaffolding |
| `--medaka_model` | `r1041_e82_400bps_hac_v4.3.0` | Medaka model; must match the basecaller model |

### Characterisation

| Parameter | Default | Description |
|---|---|---|
| `--phabox2_db` | *(required)* | PhaBOX2 database directory |
| `--phabox2_min_len` | `1000` | Minimum contig length for PhaBOX2 |
| `--phabox2_skip_phamer` | `Y` | Skip PhaMer virus identification |
| `--pharokka_db` | *(required)* | Pharokka database directory (version-matched) |

### Host branch

| Parameter | Default | Description |
|---|---|---|
| `--assemble_host` | `false` | Enable host recovery and typing |
| `--pool_host` | `false` | Additionally assemble pooled host reads of all samples |
| `--host_genome_size` | `5.5m` | Expected host genome size |
| `--host_min_reads` | `1000` | Skip host assembly below this number of host reads |
| `--host_min_depth` | `10` | Flag samples below this expected host depth as `LOW` |
| `--kaptive_k_db` / `--kaptive_o_db` | `kpsc_k` / `kpsc_o` | Kaptive databases |
| `--mlst_scheme` | `klebsiella` | PubMLST scheme; `false` for automatic detection |

### Resources

| Parameter | Default | Description |
|---|---|---|
| `--threads` | `8` | CPUs requested per multithreaded process |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `16.GB` / `48.h` | Upper limits applied to every process; raise on larger machines or HPC |

---

## 8. Output

```
results/
├── 01_porechop_abi/
│   ├── <sample>_trimmed.fastq.gz          adapter-trimmed reads
│   └── <sample>_porechop.log
├── 02_filtlong/
│   ├── <sample>_filtered.fastq.gz         reads used for assembly and polishing
│   └── <sample>_filtlong.log              read counts before/after filtering
├── 03_flye_assembly/<sample>/
│   ├── assembly.fasta                     draft assembly
│   ├── assembly_graph.gfa                 assembly graph (inspect in Bandage)
│   └── assembly_info.txt                  per-contig length, coverage, circularity, repeat, multiplicity
├── 04_medaka/<sample>/
│   └── <sample>_polished.fasta            polished assembly
├── 05_dnaapler/
│   └── <sample>_reoriented.fasta          FINAL PHAGE GENOME (starts at terL)
├── 06_assembly_stats/
│   └── <sample>_assembly_stats.tsv        contigs, length, largest contig, N50, N90, GC %
├── 07_phabox2/<sample>/<sample>_phabox2/
│   └── final_prediction/                  lineage, PhaGCN score, genus cluster, lifestyle, host and scores
├── 08_pharokka/<sample>/pharokka_out/
│   ├── <sample>.gbk                       annotated genome (GenBank)
│   ├── <sample>.gff                       annotation (GFF3)
│   ├── <sample>.tbl                       NCBI feature table (for submission)
│   ├── <sample>_cds_final_merged_output.tsv   per-gene coordinates, strand, PHROG, function, VFDB/CARD hits
│   ├── <sample>_cds_functions.tsv         gene counts per PHROG category + tRNA/CRISPR/tmRNA/VFDB/CARD
│   ├── <sample>_length_gc_cds_density.tsv length, GC, coding density
│   ├── <sample>_top_hits_mash_inphared.tsv    closest INPHARED genome, Mash distance, taxonomy, isolation host
│   ├── phanotate.faa / phanotate.ffn      predicted proteins / genes
│   └── terL.faa / terL.ffn                large terminase subunit (phylogenetic marker)
├── 09_host/                               (--assemble_host)
│   ├── <sample>/
│   │   ├── <sample>_host_reads.fastq.gz   non-phage reads
│   │   ├── <sample>_host_reads.stats.tsv  host reads, host bases, expected depth
│   │   ├── flye/<sample>_host_assembly.fasta
│   │   ├── <sample>_host_polished.fasta   polished host assembly
│   │   ├── kaptive/                       K- and O-locus reports (+ locus plots)
│   │   └── <sample>_mlst.tsv
│   ├── pooled/                            (--pool_host)
│   ├── host_typing_summary.tsv            per-sample assembly metrics, KL/OL locus, confidence, ST, alleles, depth, flag
│   └── host_concordance.txt               samples grouped by host profile
├── summary/assembly_summary.tsv           phage assembly statistics for all samples
└── pipeline_report.html, pipeline_timeline.html, pipeline_dag.html, pipeline_trace.txt
```

---

## 9. Interpreting the results

### Genome completeness

A complete phage genome is typically recovered as **one contig of the expected length with high, uniform coverage** (`assembly_info.txt`). Circularity should not be used as the criterion on its own:

* Many phages have **linear** genomes with direct terminal repeats (DTRs) or cohesive ends; Flye may collapse a DTR and report the contig as circular (`circ. = Y`), or report a linear genome correctly as `N`.
* Confirm completeness by comparison with the closest reference genome (Pharokka Mash hit) and, where relevant, by inspecting read pile-ups at the contig ends for terminal repeats.

### Taxonomy and novelty

* Pharokka `top_hits_mash_inphared.tsv`: a Mash distance close to 0 with ~1,000/1,000 matching hashes indicates a near-identical genome already in GenBank; larger distances (or no hit at `--mash_distance 0.2`) indicate a more divergent or novel phage.
* PhaBOX2 `GenusCluster = known_genus` places the phage in an existing ICTV genus. Family names may differ between PhaBOX2 and INPHARED because they follow different ICTV releases.

### Therapeutic suitability

PhaTYP `virulent`, absence of integrase/excisionase genes (Pharokka category *integration and excision* = 0) and absence of VFDB/CARD hits are consistent with a strictly lytic phage free of known virulence and resistance genes.

### Host typing

Use `depth_note` in `host_typing_summary.tsv`: `LOW` samples support species identification only; rely on K/O and ST calls from `OK` samples or from the pooled assembly.

---

## 10. Validation

The workflow was run end to end with `-profile conda` (environments built from bioconda) on a public dataset not used during development: *Klebsiella* phage **Mera**, BioProject PRJNA914245, ONT run **SRR22972385** (MinION, R9.4.1, PCR-based library).

| Metric | Published genome (OQ579025) | phage-assembly-nf |
|---|---|---|
| Contigs | 1 | 1 |
| Length | 41,400 bp | 41,397 bp |
| Coverage | — | 373× (after Filtlong) |
| INPHARED Mash top hit | — | OQ579025, 997/1,000 hashes, distance 7.2 × 10⁻⁵ |
| Genus (PhaGCN) | *Przondovirus* | *Przondovirus* (score 1.00) |
| Lifestyle (PhaTYP) | — | virulent (1.0) |
| Host genus (CHERRY) | *Klebsiella* | *Klebsiella* (0.84) |
| CDS / coding density | 51 / 92.5 % | 58 (PHANOTATE) / 94.4 % |

These data were generated on R9.4.1 and polished with the default R10.4.1 Medaka model, so the 3 bp length difference is an upper bound on the error expected with a matching model. The host branch recovered 9,587 non-phage reads (~2.5× of a 5.5 Mb genome), producing a fragmented 693 kb assembly (61 contigs, 58.4 % GC) that was flagged `LOW`; BLAST of host contigs identified *Klebsiella*.

---

## 11. Computational requirements

* **Phage branch:** typically < 8 GB RAM and minutes to an hour per isolate on a laptop, depending on depth (use `--target_bases` for very deep samples).
* **Host branch:** Flye on a ~5 Mb bacterial genome requests up to 16 GB RAM.
* **First run with `-profile conda`:** environment creation can take 30–60 min; environments are cached in `work/conda/` and reused.
* **Disk:** Pharokka database ~2 GB; PhaBOX2 database several GB.

---

## 12. Troubleshooting

| Message | Cause | Fix |
|---|---|---|
| `Cannot find Java or it's a wrong version` | Java < 17 in the active environment | Use the dedicated Nextflow environment (section 4.1) |
| `Process requirement exceeds available memory` | Machine has less memory than a step requests | `--max_memory 8.GB` (or the available amount) |
| `PHROGs Database file VERSION_x_x_x is missing` | Pharokka database from a different release | Install the database with Pharokka 1.10.1 (section 4.4) |
| `--phabox2_db directory not found` / `is empty` | Wrong database path | Correct the path; checked at launch |
| `TypeError: Invalid value 'known_genus' for dtype 'int64'` | PhaBOX2 run with pandas 3 | Use the pipeline's environment (pandas 2.2 pinned) |
| `Unable to acquire lock on session` | Another Nextflow run is active in the same directory, or a stale lock | Stop the other run, or remove the reported `LOCK` file if no process holds it |
| Host branch flagged `LOW` | Little host DNA in the lysate (DNase-treated) | Expected; use for species ID only, or `--pool_host` if isolates share a host |

---

## 13. Software versions

Pinned in `nextflow.config` for the `conda` and `docker` profiles.

| Tool | Version | Tool | Version |
|---|---|---|---|
| Porechop_ABI | 0.5.0 | Pharokka | 1.10.1 |
| Filtlong | 0.3.1 | minimap2 | 2.28 |
| Flye | 2.9.6 | SAMtools | 1.20 |
| Medaka | 2.1.1 | Kaptive | 3.1.0 |
| Dnaapler | 1.4.0 | mlst | 2.33.1 |
| PhaBOX2 | 2.1.13 (pandas 2.2) | Nextflow | ≥ 24.04 (tested 25.10.0) |

Tested on macOS (Intel) with Nextflow 25.10.0 and OpenJDK 17.

---

## 14. Citation

If you use this workflow, please cite:

> Algarni, A. T. (2026). *phage-assembly-nf: Nextflow pipeline for assembly, polishing and annotation of bacteriophage genomes and typing of their propagation host from Oxford Nanopore reads*. Zenodo. https://doi.org/10.5281/zenodo.23300714

This concept DOI always resolves to the latest version. To cite a specific release, use its version DOI: v1.0.0 — [10.5281/zenodo.23300715](https://doi.org/10.5281/zenodo.23300715); v1.0.1 — [10.5281/zenodo.23301839](https://doi.org/10.5281/zenodo.23301839).

Please also cite the tools and databases the workflow depends on (section 15).

## 15. References

**Workflow management**

1. Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nat Biotechnol.* 2017;35:316–319. https://doi.org/10.1038/nbt.3820

**Read processing, assembly and polishing**

2. Bonenfant Q, Noé L, Touzet H. Porechop_ABI: discovering unknown adapters in Oxford Nanopore Technology sequencing reads for downstream trimming. *Bioinform Adv.* 2023;3(1):vbac085. https://doi.org/10.1093/bioadv/vbac085
3. Wick RR. Filtlong. https://github.com/rrwick/Filtlong
4. Kolmogorov M, Yuan J, Lin Y, Pevzner PA. Assembly of long, error-prone reads using repeat graphs. *Nat Biotechnol.* 2019;37:540–546. https://doi.org/10.1038/s41587-019-0072-8
5. Oxford Nanopore Technologies. Medaka. https://github.com/nanoporetech/medaka
6. Bouras G, Grigson SR, Papudeshi B, Mallawaarachchi V, Roach MJ. Dnaapler: a tool to reorient circular microbial genomes. *J Open Source Softw.* 2024;9(93):5968. https://doi.org/10.21105/joss.05968

**Phage characterisation**

7. Shang J, Peng C, Liao H, Tang X, Sun Y. PhaBOX: a web server for identifying and characterizing phage contigs in metagenomic data. *Bioinform Adv.* 2023;3(1):vbad101. https://doi.org/10.1093/bioadv/vbad101
8. Shang J, Jiang J, Sun Y. Bacteriophage classification for assembled contigs using graph convolutional network. *Bioinformatics.* 2021;37(Suppl 1):i25–i33. https://doi.org/10.1093/bioinformatics/btab293
9. Shang J, Tang X, Sun Y. PhaTYP: predicting the lifestyle for bacteriophages using BERT. *Brief Bioinform.* 2023;24(1):bbac487. https://doi.org/10.1093/bib/bbac487
10. Shang J, Sun Y. CHERRY: a Computational metHod for accuratE pRediction of virus–pRokarYotic interactions using a graph encoder–decoder model. *Brief Bioinform.* 2022;23(5):bbac182. https://doi.org/10.1093/bib/bbac182
11. Bouras G, Nepal R, Houtak G, Psaltis AJ, Wormald PJ, Vreugde S. Pharokka: a fast scalable bacteriophage annotation tool. *Bioinformatics.* 2023;39(1):btac776. https://doi.org/10.1093/bioinformatics/btac776
12. McNair K, Zhou C, Dinsdale EA, Souza B, Edwards RA. PHANOTATE: a novel approach to gene identification in phage genomes. *Bioinformatics.* 2019;35(22):4537–4542. https://doi.org/10.1093/bioinformatics/btz265
13. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol.* 2017;35:1026–1028. https://doi.org/10.1038/nbt.3988
14. Larralde M, Zeller G. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics.* 2023;39(5):btad214. https://doi.org/10.1093/bioinformatics/btad214
15. Terzian P, Olo Ndela E, Galiez C, Lossouarn J, Pérez Bucio RE, Mom R, et al. PHROG: families of prokaryotic virus proteins clustered using remote homology. *NAR Genom Bioinform.* 2021;3(3):lqab067. https://doi.org/10.1093/nargab/lqab067
16. Liu B, Zheng D, Zhou S, Chen L, Yang J. VFDB 2022: a general classification scheme for bacterial virulence factors. *Nucleic Acids Res.* 2022;50(D1):D912–D917. https://doi.org/10.1093/nar/gkab1107
17. Alcock BP, Huynh W, Chalil R, Smith KW, Raphenya AR, Wlodarski MA, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* 2023;51(D1):D690–D699. https://doi.org/10.1093/nar/gkac920
18. Chan PP, Lin BY, Mak AJ, Lowe TM. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res.* 2021;49(16):9077–9096. https://doi.org/10.1093/nar/gkab688
19. Laslett D, Canback B. ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences. *Nucleic Acids Res.* 2004;32(1):11–16. https://doi.org/10.1093/nar/gkh152
20. Bland C, Ramsey TL, Sabree F, Lowe M, Brown K, Kyrpides NC, Hugenholtz P. CRISPR Recognition Tool (CRT): a tool for automatic detection of clustered regularly interspaced palindromic repeats. *BMC Bioinformatics.* 2007;8:209. https://doi.org/10.1186/1471-2105-8-209 (MinCED: https://github.com/ctSkennerton/minced)
21. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 2016;17:132. https://doi.org/10.1186/s13059-016-0997-x
22. Cook R, Brown N, Redgwell T, Rihtman B, Barnes M, Clokie M, et al. INfrastructure for a PHAge REference Database: identification of large-scale biases in the current collection of cultured phage genomes. *PHAGE.* 2021;2(4):214–223. https://doi.org/10.1089/phage.2021.0007

**Host branch**

23. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics.* 2018;34(18):3094–3100. https://doi.org/10.1093/bioinformatics/bty191
24. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. https://doi.org/10.1093/gigascience/giab008
25. Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and accurate in silico antigen typing with Kaptive 3. *Microb Genom.* 2025;11(6):001428. https://doi.org/10.1099/mgen.0.001428
26. Seemann T. mlst. https://github.com/tseemann/mlst
27. Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. *Wellcome Open Res.* 2018;3:124. https://doi.org/10.12688/wellcomeopenres.14826.1

**Software distribution**

28. Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nat Methods.* 2018;15:475–476. https://doi.org/10.1038/s41592-018-0046-7

## Author

Abdulmalek T. Algarni ([ORCID 0000-0002-2872-3724](https://orcid.org/0000-0002-2872-3724)) — King Abdulaziz City for Science and Technology (KACST)

## License

MIT — see `LICENSE`.
