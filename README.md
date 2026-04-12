# BSA-seq for Arabidopsis (Mutant Mapping)

A BSA-seq (Bulked Segregant Analysis) fine-mapping pipeline for *Arabidopsis thaliana* **mutants**.
This pipeline compares a **mutant pool** (mut) against a **wild-type pool** (wt) to identify candidate causal mutations.

Designed for users without a bioinformatics background — one command gets you from raw FASTQ to a candidate mutation list.

---

## How It Works

In a BSA-seq experiment:

- **Mutant pool (mut)** — plants with the mutant phenotype (enriched for the causal allele)
- **Wild-type pool (wt)** — plants with wild-type phenotype (wild-type allele)

The pipeline calculates a **Δ(delta)-SNP-index** at each variant site:

```
ΔSNP-index = SNP-index(mut pool) − SNP-index(wt pool)
```

Sites with a high Δ SNP-index (close to 1.0) in the candidate region are likely linked to the causal mutation.

---

## Pipeline Overview

```
Raw FASTQ (mut pool + wt pool)
    │
    ├─ Step 1: Quality control (Trimmomatic)
    ├─ Step 2: Alignment to TAIR10 (BWA-MEM)
    ├─ Step 3: MarkDuplicates + BQSR (GATK)
    ├─ Step 4: Variant calling (GATK HaplotypeCaller)
    ├─ Step 5: Variant filtration (hard filtering)
    ├─ Step 6: Δ SNP-index calculation → candidate site selection
    ├─ Step 7: Functional annotation (snpEff)
    └─ Step 8: Mutation classification output
```

---

## Before You Start

Two things must be prepared **manually** before running any script:

### 1. Organise your FASTQ files

```
RawData/
├── mut_pool/       ← mutant phenotype sequencing data
│   ├── sample_1.fq.gz
│   └── sample_2.fq.gz
└── wt_pool/        ← wild-type sequencing data
    ├── sample_1.fq.gz
    └── sample_2.fq.gz
```

> File names must end with `_1.fq.gz` (R1) and `_2.fq.gz` (R2).

### 2. Download the reference genome (one-time, ~120 MB)

```bash
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

Place the `.fa` file in the same directory as the scripts.

---

## Quick Start

### First-time use — auto-install environment

```bash
bash setup_and_run.sh
```

The script will:
1. Automatically install conda/micromamba (if not present)
2. Create the `BSA_seq` environment and install all required tools
3. Download the snpEff Arabidopsis database
4. Prompt you to **select your mutant and wild-type pool folders**
5. Prompt you to **enter the candidate chromosomal interval** (e.g. Chr1: 10,000,000 – 15,000,000)
6. Launch the full pipeline in the background

> The interactive prompts take less than 1 minute. The analysis itself takes ~2–4 hours.

### Subsequent runs — environment already set up

```bash
bash quick_start.sh
```

---

## Interactive Prompts Explained

When `setup_and_run.sh` reaches the data setup step, it will show:

```
[INFO] Select the mutant pool folder:
1) mut_pool
2) wt_pool
3) Enter manually
#?          ← type 1 and press Enter

[INFO] Select the wild-type pool folder:
1) mut_pool
2) wt_pool
3) Enter manually
#?          ← type 2 and press Enter

Chromosome (default: 1):            ← press Enter to use default, or type a number
Start position (default: 10000000): ← type your region start
End position   (default: 15000000): ← type your region end
```

> **Important:** Mutant pool must always be selected first. The Δ SNP-index direction depends on sample order.

---

## Output Files

```
bsa_output/result/
├── candidate_mutations.csv        ⭐ Key output: candidate sites filtered by Δ SNP-index
├── all_variants_ratio.csv         All variant sites with SNP-index values
├── snpEff_stats.html              Annotation summary (open in browser)
├── fixed_missense_mutations.csv   Homozygous missense mutations in both pools
├── all_missense_mutations.csv
├── fixed_frameshift_mutations.csv
├── all_frameshift_mutations.csv
├── fixed_nonsense_mutations.csv
└── all_nonsense_mutations.csv
```

**`fixed_*` files** contain mutations that are homozygous alt (1/1) in all samples — the most likely causal candidates.

---

## File Descriptions

| File | Description |
|------|-------------|
| `setup_and_run.sh` | One-click install + run (first-time use) |
| `quick_start.sh` | Quick launch (environment already ready) |
| `run_pipeline.sh` | Main analysis pipeline (called by the above scripts) |
| `config.sh` | Configuration file (auto-generated; can also be edited manually) |
| `calc_ratio.py` | Δ SNP-index calculation and candidate site filtering |
| `extract_mutations.py` | Extract and classify mutations from snpEff-annotated VCF |
| `annotate_mutations.py` | Optional: custom GFF+FASTA annotation for independent validation |

---

## Important Notes

- **Sample order matters:** Mutant pool must be **first**, wild-type pool **second**. This determines the sign of Δ SNP-index.
- **BQSR:** Requires a known-variants VCF. Leave `KNOWN_VCF` empty in `config.sh` to skip (minor impact on results).
- **`annotate_mutations.py`** is an optional independent validation tool. For routine analysis, `extract_mutations.py` is sufficient.
- **Don't know your target interval?** Set a broad range first (e.g. start=1, end=30000000), then narrow down based on the Δ SNP-index distribution in `all_variants_ratio.csv`.
