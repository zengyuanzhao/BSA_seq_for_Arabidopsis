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

Sites with a high positive Δ SNP-index (close to 1.0) in the candidate region are enriched in the mutant pool and are often the first candidates to inspect. Negative Δ SNP-index sites are reported separately because they can be useful for quality control or alternative segregation patterns.

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
    ├─ Step 7: Whole-genome functional annotation (snpEff)
    └─ Step 8: Mutation classification output
```

---

## Before You Start

Only one thing must be prepared **manually** before running any script:

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

### 2. Reference genome (auto-downloaded)

The `setup_and_run.sh` script will automatically create a `reference/` folder and download the TAIR10 genome if it is not already present. You can also place the file manually at:

```
reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```

---

## Quick Start

### First-time use — auto-install environment

```bash
bash setup_and_run.sh
```

The script will:
1. Automatically install conda/micromamba (if not present)
2. Create the `BSA_seq` environment and install all required tools
3. Download the snpEff `Arabidopsis_thaliana` database
4. Create a `reference/` folder and download the TAIR10 genome (if missing)
5. Prompt you to **select your mutant and wild-type pool folders**
6. Prompt you to **enter the candidate chromosomal interval** (e.g. Chr1: 10,000,000 – 15,000,000)
7. Launch the full pipeline in the background

> The interactive prompts take less than 1 minute. The analysis itself takes ~2–4 hours.

### Subsequent runs — environment already set up

```bash
bash quick_start.sh
```

`quick_start.sh` validates the required tools before launching the long-running job, so missing commands such as `snpEff`, `tabix`, or `bgzip` fail early.

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
├── candidate_positive_mutations.csv Mutant-pool-enriched candidate sites
├── candidate_negative_mutations.csv WT-pool-enriched candidate sites
├── all_genome_mutations.csv      Whole-genome compact table:
│                                  CHROM/POS/REF/ALT/QUAL/read counts/frequencies/type/GENE_ID
├── all_annotated_mutations.csv    All snpEff annotation types with gene IDs
├── candidate_annotated_mutations.csv Candidate sites with all snpEff annotation types
├── candidate_functional_mutations.csv Candidate sites with missense/frameshift/nonsense effects
├── all_variants_ratio.csv         All variant sites with SNP-index values
├── snpEff_stats.html              Annotation summary (open in browser)
├── fixed_missense_mutations.csv   Homozygous missense mutations in both pools
├── all_missense_mutations.csv
├── fixed_frameshift_mutations.csv
├── all_frameshift_mutations.csv
├── fixed_nonsense_mutations.csv
└── all_nonsense_mutations.csv
```

**`all_genome_mutations.csv`** is the compact whole-genome table with columns `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, `mut_ref`, `mut_alt`, `mut_freq`, `wt_ref`, `wt_alt`, `wt_freq`, `mut-wt_freq`, `variant`, and `GENE_ID`. `mut_freq` and `wt_freq` are calculated from the selected ALT allele depth as `alt / (ref + alt)`, and `mut-wt_freq = mut_freq - wt_freq`.

**`all_annotated_mutations.csv`** contains every snpEff annotation row from the annotated VCF, including mutation type (`MUTATION_TYPE`) and `GENE_ID`. A single `CHROM/POS/REF/ALT` may appear multiple times because snpEff can annotate one variant against multiple genes, transcripts, or effect types; this is intentional so gene IDs are not lost.

**`candidate_positive_mutations.csv`** is usually the first table to inspect for a mutant-pool-enriched causal allele. **`candidate_annotated_mutations.csv`** intersects Δ SNP-index candidates with all snpEff annotation types, while **`candidate_functional_mutations.csv`** keeps only missense / frameshift / nonsense candidate effects.

**`fixed_*` files** contain variants that are homozygous alt (1/1) in all samples. These are useful for tracking shared background variants, but they are not automatically the most likely causal mutations.

---

## File Descriptions

| File | Description |
|------|-------------|
| `setup_and_run.sh` | One-click install + run (first-time use) |
| `quick_start.sh` | Quick launch (environment already ready; validates required commands before starting) |
| `run_pipeline.sh` | Main analysis pipeline; Step 7 performs whole-genome snpEff annotation on `final_variants.vcf.gz` |
| `config.sh` | Configuration file (auto-generated; can also be edited manually) |
| `calc_ratio.py` | Δ SNP-index calculation and candidate site filtering; writes all / positive / negative candidate tables |
| `extract_mutations.py` | Extract mutation tables from whole-genome snpEff-annotated VCF, including `all_genome_mutations.csv` |
| `annotate_mutations.py` | Optional: custom GFF+FASTA annotation for independent validation |
| `.gitignore` | Excludes large data and output files from version control |

---

## Important Notes

- **Sample order matters:** Mutant pool must be **first**, wild-type pool **second**. This determines the sign of Δ SNP-index.
- **Candidate direction:** `candidate_mutations.csv` keeps both positive and negative Δ sites by default. For typical mutant-pool enrichment, inspect `candidate_positive_mutations.csv` first.
- **Whole-genome annotation:** Step 7 annotates `final_variants.vcf.gz` directly, not only the fine-mapping interval. This is required for `all_genome_mutations.csv`.
- **snpEff database:** The scripts prefer the exact `Arabidopsis_thaliana` database. Avoid using test databases such as `testAthalianaTair10`.
- **BQSR:** Requires a known-variants VCF. Leave `KNOWN_VCF` empty in `config.sh` to skip (minor impact on results).
- **`annotate_mutations.py`** is an optional independent validation tool. For routine analysis, `extract_mutations.py` is sufficient.
- **Don't know your target interval?** Set a broad range first (e.g. start=1, end=30000000), then narrow down based on the Δ SNP-index distribution in `all_variants_ratio.csv`.
- **Large files are ignored by Git:** `RawData/`, `reference/`, `bsa_output*/` and `run*.log` are listed in `.gitignore` and will not be pushed to the repository.

---

## Lightweight Validation

After editing scripts, run:

```bash
bash -n setup_and_run.sh quick_start.sh run_pipeline.sh
python3 -m py_compile calc_ratio.py extract_mutations.py annotate_mutations.py
git diff --check
```
