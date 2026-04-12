#!/bin/bash
# run_pipeline.sh - Main BSA-seq analysis pipeline
# Do NOT run this script directly.
# Use quick_start.sh (environment already set up) or setup_and_run.sh (first time).

set -euo pipefail

# Load configuration (shell variable format)
if [ ! -f "config.sh" ]; then
    echo "[ERROR] config.sh not found. Please run setup_and_run.sh first."
    exit 1
fi
source config.sh

# ---- Validate required config variables ----
for var in MUT_R1 MUT_R2 WT_R1 WT_R2; do
    if [[ -z "${!var}" ]]; then
        echo "[ERROR] ${var} is not set in config.sh. Please fill in your FASTQ file paths."
        exit 1
    fi
done

mkdir -p "${OUTDIR}/qc" "${OUTDIR}/bam" "${OUTDIR}/vcf" "${OUTDIR}/result"
LOG="${OUTDIR}/pipeline.log"
exec > >(tee -a "${LOG}") 2>&1

echo "====== [$(date)] Pipeline started ======"
echo "  Mutant pool: ${MUT_R1}"
echo "  WT pool:     ${WT_R1}"
echo "  Target region:   ${FINE_CHR}:${FINE_START}-${FINE_END}"
echo ""

# ---- Check and build reference genome indices ----
echo "[Prep] Checking reference genome indices..."
REF_DICT="${REF%.fa}.dict"
[[ "${REF}" == *.fasta ]] && REF_DICT="${REF%.fasta}.dict"

if [[ ! -f "${REF}.fai" ]]; then
    echo "  Creating samtools index (.fai)..."
    samtools faidx "${REF}"
    echo "  Done: .fai"
else
    echo "  OK: .fai already exists"
fi

if [[ ! -f "${REF_DICT}" ]]; then
    echo "  Creating GATK sequence dictionary (.dict)..."
    gatk CreateSequenceDictionary -R "${REF}" -O "${REF_DICT}"
    echo "  Done: .dict"
else
    echo "  OK: .dict already exists"
fi

# ---- Step 1: Quality control ----
echo "[Step 1] Quality control with Trimmomatic..."
for SAMPLE in mut wt; do
    R1="${OUTDIR}/qc/${SAMPLE}_1.fq.gz"
    R2="${OUTDIR}/qc/${SAMPLE}_2.fq.gz"
    UNPAIRED_R1="${OUTDIR}/qc/${SAMPLE}_1_unpaired.fq.gz"
    UNPAIRED_R2="${OUTDIR}/qc/${SAMPLE}_2_unpaired.fq.gz"
    if [[ ${SAMPLE} == "mut" ]]; then
        IN1=${MUT_R1}; IN2=${MUT_R2}
    else
        IN1=${WT_R1};  IN2=${WT_R2}
    fi
    trimmomatic PE -threads ${THREADS} ${IN1} ${IN2} \
        ${R1} ${UNPAIRED_R1} ${R2} ${UNPAIRED_R2} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36
done

# ---- Step 2: Alignment ----
echo "[Step 2] BWA-MEM alignment..."
if [[ ! -f "${REF}.bwt" ]]; then
    echo "  Building BWA index (this may take several minutes)..."
    bwa index "${REF}"
    echo "  Done: BWA index"
else
    echo "  OK: BWA index already exists"
fi
for SAMPLE in mut wt; do
    R1="${OUTDIR}/qc/${SAMPLE}_1.fq.gz"
    R2="${OUTDIR}/qc/${SAMPLE}_2.fq.gz"
    bwa mem -t ${THREADS} \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "${REF}" ${R1} ${R2} | \
        samtools sort -@ ${THREADS} -o "${OUTDIR}/bam/${SAMPLE}_sorted.bam"
    samtools index "${OUTDIR}/bam/${SAMPLE}_sorted.bam"
done

# ---- Step 3: Mark duplicates + BQSR ----
echo "[Step 3] MarkDuplicates + BQSR..."
for SAMPLE in mut wt; do
    BAM="${OUTDIR}/bam/${SAMPLE}_sorted.bam"
    MKDUP="${OUTDIR}/bam/${SAMPLE}_markdup.bam"
    RECAL="${OUTDIR}/bam/${SAMPLE}_recal.bam"

    gatk MarkDuplicates -I "${BAM}" -O "${MKDUP}" \
        -M "${OUTDIR}/bam/${SAMPLE}_metrics.txt"
    samtools index "${MKDUP}"

    if [[ -n "${KNOWN_VCF}" && -f "${KNOWN_VCF}" ]]; then
        gatk BaseRecalibrator -I "${MKDUP}" -R "${REF}" \
            --known-sites "${KNOWN_VCF}" \
            -O "${OUTDIR}/bam/${SAMPLE}_recal.table"
        gatk ApplyBQSR -I "${MKDUP}" -R "${REF}" \
            --bqsr-recal-file "${OUTDIR}/bam/${SAMPLE}_recal.table" \
            -O "${RECAL}"
        samtools index "${RECAL}"
    else
        echo "  Skipping BQSR (KNOWN_VCF not set)"
        (cd "${OUTDIR}/bam" && ln -sf "$(basename ${MKDUP})" "$(basename ${RECAL})")
        samtools index "${RECAL}"
    fi
done

# ---- Step 4: HaplotypeCaller ----
echo "[Step 4] Calling SNPs/Indels (HaplotypeCaller)..."
for SAMPLE in mut wt; do
    gatk HaplotypeCaller -R "${REF}" \
        -I "${OUTDIR}/bam/${SAMPLE}_recal.bam" \
        -O "${OUTDIR}/vcf/${SAMPLE}.g.vcf.gz" \
        -ERC GVCF --sample-name "${SAMPLE}"
done

gatk CombineGVCFs -R "${REF}" \
    -V "${OUTDIR}/vcf/mut.g.vcf.gz" \
    -V "${OUTDIR}/vcf/wt.g.vcf.gz" \
    -O "${OUTDIR}/vcf/combined.g.vcf.gz"

gatk GenotypeGVCFs -R "${REF}" \
    -V "${OUTDIR}/vcf/combined.g.vcf.gz" \
    -O "${OUTDIR}/vcf/raw_variants.vcf.gz"

# ---- Step 5: Variant filtration (hard filtering) ----
echo "[Step 5] Variant filtration..."
gatk SelectVariants -V "${OUTDIR}/vcf/raw_variants.vcf.gz" \
    -R "${REF}" --select-type-to-include SNP \
    -O "${OUTDIR}/vcf/raw_snps.vcf.gz"
gatk VariantFiltration -R "${REF}" -V "${OUTDIR}/vcf/raw_snps.vcf.gz" \
    --filter-expression "QD < 2.0"  --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    -O "${OUTDIR}/vcf/filtered_snps.vcf.gz"

gatk SelectVariants -V "${OUTDIR}/vcf/raw_variants.vcf.gz" \
    -R "${REF}" --select-type-to-include INDEL \
    -O "${OUTDIR}/vcf/raw_indels.vcf.gz"
gatk VariantFiltration -R "${REF}" -V "${OUTDIR}/vcf/raw_indels.vcf.gz" \
    --filter-expression "QD < 2.0"   --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    -O "${OUTDIR}/vcf/filtered_indels.vcf.gz"

gatk MergeVcfs \
    -I "${OUTDIR}/vcf/filtered_snps.vcf.gz" \
    -I "${OUTDIR}/vcf/filtered_indels.vcf.gz" \
    -O "${OUTDIR}/vcf/final_variants.vcf.gz"

# ---- Step 6: Delta-SNP-index calculation ----
echo "[Step 6] Calculating delta-SNP-index and selecting candidate variants..."
python3 calc_ratio.py \
    --vcf   "${OUTDIR}/vcf/final_variants.vcf.gz" \
    --chr   "${FINE_CHR}" \
    --start "${FINE_START}" \
    --end   "${FINE_END}" \
    --outdir "${OUTDIR}/result"

# ---- Step 7: snpEff annotation (fine-mapping region only) ----
echo "[Step 7] snpEff variant annotation (region: ${FINE_CHR}:${FINE_START}-${FINE_END})..."
bcftools view -r "${FINE_CHR}:${FINE_START}-${FINE_END}" \
    "${OUTDIR}/vcf/final_variants.vcf.gz" \
    -O z -o "${OUTDIR}/vcf/region_variants.vcf.gz"
tabix "${OUTDIR}/vcf/region_variants.vcf.gz"

_DB_LIST=$(snpEff databases 2>/dev/null || true)
SNPEFF_DB=$(echo "${_DB_LIST}" | grep -i 'TAIR10' | head -1 | awk '{print $1}' || true)
if [ -z "${SNPEFF_DB}" ]; then
    echo "  TAIR10 database not found locally, attempting download..."
    if ! snpEff download TAIR10.31; then
        echo "[ERROR] Failed to download TAIR10.31 database."
        echo "        Please check your network connection, or manually run:"
        echo "            snpEff download TAIR10.31"
        exit 1
    fi
    SNPEFF_DB="TAIR10.31"
fi
echo "  Using database: ${SNPEFF_DB}"

snpEff ann -v \
    -stats "${OUTDIR}/result/snpEff_stats.html" \
    "${SNPEFF_DB}" \
    "${OUTDIR}/vcf/region_variants.vcf.gz" \
    > "${OUTDIR}/vcf/annotated_variants.vcf" 2>> "${LOG}"

bgzip -f "${OUTDIR}/vcf/annotated_variants.vcf"
tabix -f "${OUTDIR}/vcf/annotated_variants.vcf.gz"
echo "  Done: snpEff annotation complete"

# ---- Step 8: Mutation classification ----
echo "[Step 8] Extracting and classifying mutations..."
python3 extract_mutations.py \
    "${OUTDIR}/vcf/annotated_variants.vcf.gz" \
    "${OUTDIR}/result"

echo ""
echo "====== [$(date)] Pipeline complete! Results in ${OUTDIR}/result/ ======"

wait
