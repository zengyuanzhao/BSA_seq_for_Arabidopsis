#!/usr/bin/env python3
"""
extract_mutations.py - Extract and classify mutations from a snpEff-annotated VCF.

Dependencies: Python 3 standard library only (no extra packages required).
Usage: python3 extract_mutations.py <annotated.vcf.gz> <output_dir> [candidate_mutations.csv]

Outputs include:
  all_genome_mutations.csv         Compact whole-genome table with allele counts and gene IDs
  all_annotated_mutations.csv      All snpEff annotations with mutation type and gene ID
  candidate_annotated_mutations.csv Candidate rows intersected with all snpEff annotation types
  candidate_functional_mutations.csv Candidate rows limited to missense/frameshift/nonsense effects
"""

import gzip
import re
import csv
import os
import sys


def load_candidate_index(candidate_file):
    """Load candidate_mutations.csv keyed by CHROM, POS, REF, ALT."""
    if not candidate_file:
        return {}
    if not os.path.isfile(candidate_file):
        print(f"WARNING: Candidate file not found, skipping candidate/effect merge: {candidate_file}")
        return {}

    candidate_index = {}
    with open(candidate_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (
                row.get('CHROM', ''),
                row.get('POS', ''),
                row.get('REF', ''),
                row.get('ALT', ''),
            )
            if all(key):
                candidate_index[key] = row
    return candidate_index


def ann_part(parts, index):
    """Return a snpEff ANN field by index, or an empty string if missing."""
    return parts[index] if len(parts) > index else ''


def add_candidate_metrics(row, candidate_row):
    """Copy SNP-index metrics from candidate_mutations.csv into an annotation row."""
    merged = row.copy()
    for field in ['mut_ratio', 'mut_dp', 'wt_ratio', 'wt_dp', 'delta_index']:
        merged[field] = candidate_row.get(field, '')
    return merged


def parse_ad(ad_str):
    """Parse a VCF AD field into integer allele depths."""
    if not ad_str or ad_str == '.':
        return []
    try:
        return [int(x) for x in ad_str.split(',') if x != '.']
    except ValueError:
        return []


def allele_metrics(ad_counts, alt_ad_index):
    """Return ref count, selected-alt count, and selected-alt frequency."""
    if not ad_counts:
        return '', '', ''
    ref_count = ad_counts[0] if len(ad_counts) > 0 else ''
    alt_count = ad_counts[alt_ad_index] if len(ad_counts) > alt_ad_index else ''
    if ref_count == '' or alt_count == '':
        return ref_count, alt_count, ''
    depth = ref_count + alt_count
    if depth == 0:
        return ref_count, alt_count, ''
    return ref_count, alt_count, round(alt_count / depth, 4)


def resolve_alt_allele(ann_allele, alt):
    """Return the ALT allele and its AD index for a snpEff ANN allele."""
    alts = alt.split(',')
    if ann_allele in alts:
        return ann_allele, alts.index(ann_allele) + 1
    return alts[0], 1


def main(vcf_file, result_dir, candidate_file=None):
    os.makedirs(result_dir, exist_ok=True)
    candidate_index = load_candidate_index(candidate_file)

    mutations = {
        'all_missense':     [],
        'fixed_missense':   [],
        'all_frameshift':   [],
        'fixed_frameshift': [],
        'all_nonsense':     [],
        'fixed_nonsense':   [],
    }
    genome_rows = []
    genome_row_keys = set()
    all_annotations = []
    candidate_annotations = []
    candidate_functional = []
    annotated_sites = set()
    candidate_annotated_sites = set()

    opener = gzip.open if vcf_file.endswith('.gz') else open
    samples = []

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                samples = headers[9:]
                continue

            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue

            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]

            info = cols[7]
            fmt  = cols[8].split(':')

            sample_gts = {}
            sample_ads = {}
            ad_idx = fmt.index('AD') if 'AD' in fmt else None
            if 'GT' in fmt:
                gt_idx = fmt.index('GT')
            else:
                gt_idx = None
            for i, sname in enumerate(samples):
                sdata = cols[9 + i].split(':')
                if gt_idx is not None:
                    gt = sdata[gt_idx] if gt_idx < len(sdata) else './.'
                    sample_gts[sname] = gt
                if ad_idx is not None:
                    ad_str = sdata[ad_idx] if ad_idx < len(sdata) else '.'
                    sample_ads[sname] = parse_ad(ad_str)

            is_fixed = (
                all(gt in ['1/1', '1|1'] for gt in sample_gts.values())
                and len(sample_gts) > 0
            )
            mut_gt = sample_gts.get('mut', '')
            wt_gt = sample_gts.get('wt', '')
            mut_fixed = mut_gt in ['1/1', '1|1']
            wt_fixed = wt_gt in ['1/1', '1|1']
            candidate_row = candidate_index.get((chrom, pos, ref, alt))

            ann_match = re.search(r'ANN=([^;]+)', info)
            if not ann_match:
                continue

            ann_list = ann_match.group(1).split(',')
            seen_annotations = set()
            seen_all_annotations = set()

            for ann in ann_list:
                parts = ann.split('|')
                if len(parts) < 2:
                    continue

                allele       = ann_part(parts, 0)
                annotation   = ann_part(parts, 1)
                impact       = ann_part(parts, 2)
                gene_name    = ann_part(parts, 3)
                gene_id      = ann_part(parts, 4)
                feature_type = ann_part(parts, 5)
                feature_id   = ann_part(parts, 6)
                biotype      = ann_part(parts, 7)
                rank         = ann_part(parts, 8)
                hgvs_c       = ann_part(parts, 9)
                hgvs_p       = ann_part(parts, 10)
                distance     = ann_part(parts, 14)
                messages     = ann_part(parts, 15)
                alt_for_row, alt_ad_index = resolve_alt_allele(allele, alt)
                mut_ref, mut_alt, mut_freq = allele_metrics(
                    sample_ads.get('mut', []), alt_ad_index
                )
                wt_ref, wt_alt, wt_freq = allele_metrics(
                    sample_ads.get('wt', []), alt_ad_index
                )
                delta_freq = ''
                if mut_freq != '' and wt_freq != '':
                    delta_freq = round(mut_freq - wt_freq, 4)

                var_info = {
                    'CHROM':              chrom,
                    'POS':                pos,
                    'REF':                ref,
                    'ALT':                alt,
                    'ALLELE':             allele,
                    'MUTATION_TYPE':      annotation,
                    'ANNOTATION':         annotation,
                    'IMPACT':             impact,
                    'GENE_NAME':          gene_name,
                    'GENE_ID':            gene_id,
                    'FEATURE_TYPE':       feature_type,
                    'FEATURE_ID':         feature_id,
                    'TRANSCRIPT_BIOTYPE': biotype,
                    'RANK':               rank,
                    'HGVS_C':             hgvs_c,
                    'HGVS_P':             hgvs_p,
                    'DISTANCE':           distance,
                    'ERRORS_WARNINGS':    messages,
                    'MUT_FIXED':          mut_fixed,
                    'WT_FIXED':           wt_fixed,
                    'IS_FIXED':           is_fixed,
                }
                for s, gt in sample_gts.items():
                    var_info[f'{s}_GT'] = gt

                genome_row = {
                    'CHROM':       chrom,
                    'POS':         pos,
                    'REF':         ref,
                    'ALT':         alt_for_row,
                    'QUAL':        cols[5],
                    'mut_ref':     mut_ref,
                    'mut_alt':     mut_alt,
                    'mut_freq':    mut_freq,
                    'wt_ref':      wt_ref,
                    'wt_alt':      wt_alt,
                    'wt_freq':     wt_freq,
                    'mut-wt_freq': delta_freq,
                    'variant':     annotation,
                    'GENE_ID':     gene_id,
                }
                genome_row_key = tuple(genome_row.values())
                if genome_row_key not in genome_row_keys:
                    genome_row_keys.add(genome_row_key)
                    genome_rows.append(genome_row)

                all_dedup_key = (
                    chrom, pos, ref, alt, allele, annotation, gene_id,
                    feature_id, hgvs_c, hgvs_p
                )
                if all_dedup_key not in seen_all_annotations:
                    seen_all_annotations.add(all_dedup_key)
                    all_annotations.append(var_info)
                    annotated_sites.add((chrom, pos, ref, alt))
                    if candidate_row:
                        candidate_annotations.append(add_candidate_metrics(var_info, candidate_row))
                        candidate_annotated_sites.add((chrom, pos, ref, alt))

                dedup_key = (chrom, pos, annotation, gene_id)
                if dedup_key in seen_annotations:
                    continue
                seen_annotations.add(dedup_key)

                is_functional = False
                if 'missense_variant' in annotation:
                    mutations['all_missense'].append(var_info)
                    is_functional = True
                    if is_fixed:
                        mutations['fixed_missense'].append(var_info)
                elif 'frameshift' in annotation:
                    mutations['all_frameshift'].append(var_info)
                    is_functional = True
                    if is_fixed:
                        mutations['fixed_frameshift'].append(var_info)
                elif 'stop_gained' in annotation or 'stop_lost' in annotation:
                    mutations['all_nonsense'].append(var_info)
                    is_functional = True
                    if is_fixed:
                        mutations['fixed_nonsense'].append(var_info)

                if is_functional and candidate_row:
                    candidate_functional.append(add_candidate_metrics(var_info, candidate_row))

    base_fields = [
        'CHROM', 'POS', 'REF', 'ALT', 'ALLELE', 'MUTATION_TYPE',
        'ANNOTATION', 'IMPACT', 'GENE_NAME', 'GENE_ID', 'FEATURE_TYPE',
        'FEATURE_ID', 'TRANSCRIPT_BIOTYPE', 'RANK', 'HGVS_C', 'HGVS_P',
        'DISTANCE', 'ERRORS_WARNINGS', 'MUT_FIXED', 'WT_FIXED', 'IS_FIXED'
    ]
    gt_fields = [f'{s}_GT' for s in samples]
    fieldnames = base_fields + gt_fields
    genome_fieldnames = [
        'CHROM', 'POS', 'REF', 'ALT', 'QUAL',
        'mut_ref', 'mut_alt', 'mut_freq',
        'wt_ref', 'wt_alt', 'wt_freq',
        'mut-wt_freq', 'variant', 'GENE_ID',
    ]
    candidate_fieldnames = (
        base_fields +
        ['mut_ratio', 'mut_dp', 'wt_ratio', 'wt_dp', 'delta_index'] +
        gt_fields
    )

    print("\n=== Mutation Classification Summary ===")
    genome_output = os.path.join(result_dir, 'all_genome_mutations.csv')
    with open(genome_output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=genome_fieldnames)
        writer.writeheader()
        if genome_rows:
            writer.writerows(genome_rows)
    print(f"  all_genome: {len(genome_rows)} rows -> {genome_output}")

    all_output = os.path.join(result_dir, 'all_annotated_mutations.csv')
    with open(all_output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        if all_annotations:
            writer.writerows(all_annotations)
    print(
        f"  all_annotated: {len(annotated_sites)} sites, "
        f"{len(all_annotations)} annotation rows -> {all_output}"
    )

    candidate_all_output = os.path.join(result_dir, 'candidate_annotated_mutations.csv')
    with open(candidate_all_output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=candidate_fieldnames)
        writer.writeheader()
        if candidate_annotations:
            writer.writerows(candidate_annotations)
    print(
        f"  candidate_annotated: {len(candidate_annotated_sites)} sites, "
        f"{len(candidate_annotations)} annotation rows -> {candidate_all_output}"
    )

    for category, var_list in mutations.items():
        output_file = os.path.join(result_dir, f'{category}_mutations.csv')
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            if var_list:
                writer.writerows(var_list)
        print(f"  {category}: {len(var_list)} sites -> {output_file}")

    candidate_output = os.path.join(result_dir, 'candidate_functional_mutations.csv')
    with open(candidate_output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=candidate_fieldnames)
        writer.writeheader()
        if candidate_functional:
            writer.writerows(candidate_functional)
    print(f"  candidate_functional: {len(candidate_functional)} sites -> {candidate_output}")


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        print("Usage: python3 extract_mutations.py <annotated.vcf.gz> <output_dir> [candidate_mutations.csv]")
        sys.exit(1)
    vcf_file   = sys.argv[1]
    result_dir = sys.argv[2]
    candidate_file = sys.argv[3] if len(sys.argv) == 4 else None
    main(vcf_file, result_dir, candidate_file)
