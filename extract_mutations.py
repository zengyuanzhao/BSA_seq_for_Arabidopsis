#!/usr/bin/env python3
"""
extract_mutations.py - Extract and classify mutations from a snpEff-annotated VCF.

Dependencies: Python 3 standard library only (no extra packages required).
Usage: python3 extract_mutations.py <annotated.vcf.gz> <output_dir>
"""

import gzip
import re
import csv
import os
import sys


def main(vcf_file, result_dir):
    os.makedirs(result_dir, exist_ok=True)

    mutations = {
        'all_missense':     [],
        'fixed_missense':   [],
        'all_frameshift':   [],
        'fixed_frameshift': [],
        'all_nonsense':     [],
        'fixed_nonsense':   [],
    }

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

            filter_val = cols[6] if len(cols) > 6 else 'PASS'
            if filter_val not in ('PASS', '.'):
                continue

            info = cols[7]
            fmt  = cols[8].split(':')

            sample_gts = {}
            if 'GT' in fmt:
                gt_idx = fmt.index('GT')
                for i, sname in enumerate(samples):
                    sdata = cols[9 + i].split(':')
                    gt = sdata[gt_idx] if gt_idx < len(sdata) else './.'
                    sample_gts[sname] = gt

            is_fixed = (
                all(gt in ['1/1', '1|1'] for gt in sample_gts.values())
                and len(sample_gts) > 0
            )

            ann_match = re.search(r'ANN=([^;]+)', info)
            if not ann_match:
                continue

            ann_list = ann_match.group(1).split(',')
            seen_annotations = set()

            for ann in ann_list:
                parts = ann.split('|')
                if len(parts) < 2:
                    continue

                annotation = parts[1]
                gene_id    = parts[4]  if len(parts) > 4  else ''
                feature_id = parts[6]  if len(parts) > 6  else ''
                hgvs_p     = parts[10] if len(parts) > 10 else ''

                dedup_key = (chrom, pos, annotation, gene_id)
                if dedup_key in seen_annotations:
                    continue
                seen_annotations.add(dedup_key)

                var_info = {
                    'CHROM':      chrom,
                    'POS':        pos,
                    'REF':        ref,
                    'ALT':        alt,
                    'GENE_ID':    gene_id,
                    'FEATURE_ID': feature_id,
                    'ANNOTATION': annotation,
                    'HGVS_P':     hgvs_p,
                    'IS_FIXED':   is_fixed,
                }
                for s, gt in sample_gts.items():
                    var_info[f'{s}_GT'] = gt

                if 'missense_variant' in annotation:
                    mutations['all_missense'].append(var_info)
                    if is_fixed:
                        mutations['fixed_missense'].append(var_info)
                elif 'frameshift' in annotation:
                    mutations['all_frameshift'].append(var_info)
                    if is_fixed:
                        mutations['fixed_frameshift'].append(var_info)
                elif 'stop_gained' in annotation or 'stop_lost' in annotation:
                    mutations['all_nonsense'].append(var_info)
                    if is_fixed:
                        mutations['fixed_nonsense'].append(var_info)

    base_fields = [
        'CHROM', 'POS', 'REF', 'ALT', 'GENE_ID', 'FEATURE_ID',
        'ANNOTATION', 'HGVS_P', 'IS_FIXED'
    ]
    gt_fields = [f'{s}_GT' for s in samples]
    fieldnames = base_fields + gt_fields

    print("\n=== Mutation Classification Summary ===")
    for category, var_list in mutations.items():
        output_file = os.path.join(result_dir, f'{category}_mutations.csv')
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            if var_list:
                writer.writerows(var_list)
        print(f"  {category}: {len(var_list)} sites -> {output_file}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 extract_mutations.py <annotated.vcf.gz> <output_dir>")
        sys.exit(1)
    vcf_file   = sys.argv[1]
    result_dir = sys.argv[2]
    main(vcf_file, result_dir)
