#!/usr/bin/env python3
"""
annotate_mutations.py - Custom variant effect annotation using GFF + FASTA (optional tool).

NOTE: The main pipeline already uses snpEff for annotation (Step 7).
      This script is an OPTIONAL independent validation tool.
      For routine analysis, use extract_mutations.py to parse snpEff results.

Known limitation:
  Amino acid position (aa_position) for minus-strand multi-exon genes is
  calculated within a single CDS exon and may be inaccurate. For precise
  positions use snpEff (Step 7 of the main pipeline).

Usage:
  python3 annotate_mutations.py \\
    --vcf  bsa_output/vcf/final_variants.vcf.gz \\
    --gff  Arabidopsis_thaliana.TAIR10.54.gff3 \\
    --ref  Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \\
    --outdir annotation_result
"""

import argparse
import gzip
import os
import csv
from collections import defaultdict

# Standard codon table
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in seq[::-1]])


def parse_gff(gff_file):
    """Parse a GFF file and build sorted CDS interval lists per chromosome."""
    cds_intervals = defaultdict(list)  # chrom -> [(start, end, strand, mrna_id, gene_id), ...]
    opener = gzip.open if gff_file.endswith('.gz') else open

    with opener(gff_file, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr = cols
            if feature != 'CDS':
                continue
            start, end = int(start), int(end)

            mrna_id = None
            for item in attr.split(';'):
                if item.startswith('Parent='):
                    for p in item[7:].split(','):
                        p = p.strip()
                        if p and not p.endswith('-Protein'):
                            mrna_id = p
                            break
                if mrna_id:
                    break

            if mrna_id:
                gene_id = mrna_id.split('.')[0]
                cds_intervals[chrom].append((start, end, strand, mrna_id, gene_id))

    for chrom in cds_intervals:
        cds_intervals[chrom].sort(key=lambda x: x[0])

    return cds_intervals


def find_all_cds(chrom, pos, cds_intervals):
    """
    Find all CDS intervals overlapping a given position.
    Returns a list to support multi-transcript genes.
    """
    if chrom not in cds_intervals:
        return []
    results = []
    for interval in cds_intervals[chrom]:
        start, end, strand, mrna_id, gene_id = interval
        if start <= pos <= end:
            results.append(interval)
        elif start > pos:
            break  # Intervals are sorted; no need to continue
    return results


def load_fasta_dict(fasta_file):
    """Load a FASTA file into a dict keyed by chromosome name."""
    sequences = {}
    current_chrom = None
    current_seq   = []
    opener = gzip.open if fasta_file.endswith('.gz') else open

    with opener(fasta_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chrom:
                    sequences[current_chrom] = ''.join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq   = []
            else:
                current_seq.append(line)
    if current_chrom:
        sequences[current_chrom] = ''.join(current_seq)
    return sequences


def get_codon_at_position(seq, chrom, pos, cds_start, cds_end, strand):
    """Return the codon, position-in-codon, and genomic codon start for a given site.

    NOTE: cds_pos is computed within this single CDS exon only.
    For minus-strand multi-exon genes the resulting aa_position may be
    inaccurate (known limitation; use snpEff for precise results).
    """
    if chrom not in seq:
        return None, None, None

    cds_pos = (pos - cds_start + 1) if strand == '+' else (cds_end - pos + 1)

    codon_start_in_cds = ((cds_pos - 1) // 3) * 3
    pos_in_codon       = (cds_pos - 1) % 3

    if strand == '+':
        codon_genome_start = cds_start + codon_start_in_cds
        codon = seq[chrom][codon_genome_start - 1 : codon_genome_start + 2]
    else:
        codon_genome_start = cds_end - codon_start_in_cds - 2
        codon = seq[chrom][codon_genome_start - 1 : codon_genome_start + 2]
        codon = reverse_complement(codon)
        pos_in_codon = 2 - pos_in_codon  # Position reverses on minus strand

    return codon.upper(), pos_in_codon, codon_genome_start


def predict_variant_effect(chrom, pos, ref, alt, cds_info, seq):
    """Predict variant effect: synonymous / missense / nonsense / frameshift / inframe indel."""
    cds_start, cds_end, strand, mrna_id, gene_id = cds_info
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1 and alt_len == 1:  # SNP
        codon, pos_in_codon, _ = get_codon_at_position(
            seq, chrom, pos, cds_start, cds_end, strand
        )
        if not codon or len(codon) != 3:
            return 'unknown', gene_id, mrna_id, None, None

        mut_codon = list(codon)
        mut_base  = alt if strand == '+' else reverse_complement(alt)
        mut_codon[pos_in_codon] = mut_base

        ref_aa = codon_table.get(codon.upper(), 'X')
        mut_aa = codon_table.get(''.join(mut_codon).upper(), 'X')
        aa_pos = ((pos - cds_start if strand == '+' else cds_end - pos) // 3) + 1

        if ref_aa == mut_aa:
            effect = 'synonymous'
        elif mut_aa == '*':
            effect = 'nonsense'
        else:
            effect = 'missense'
        return effect, gene_id, mrna_id, aa_pos, (ref_aa, mut_aa)

    else:  # Indel
        len_change = alt_len - ref_len
        if len_change % 3 == 0:
            effect = 'inframe_insertion' if len_change > 0 else 'inframe_deletion'
        else:
            effect = 'frameshift'
        return effect, gene_id, mrna_id, None, (None, None)


def parse_vcf(vcf_file):
    """Parse a VCF file and return a list of variant dicts plus sample names."""
    variants = []
    samples  = []
    opener   = gzip.open if vcf_file.endswith('.gz') else open

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

            chrom, pos, vid, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
            filter_val = cols[6] if len(cols) > 6 else 'PASS'
            if filter_val not in ('PASS', '.'):
                continue

            fmt = cols[8].split(':')
            sample_gts = {}
            if 'GT' in fmt:
                gt_idx = fmt.index('GT')
                for i, sname in enumerate(samples):
                    sdata = cols[9 + i].split(':')
                    gt = sdata[gt_idx] if gt_idx < len(sdata) else './.'
                    sample_gts[sname] = gt

            variants.append({
                'CHROM': chrom, 'POS': pos, 'ID': vid,
                'REF': ref, 'ALT': alt,
                'FILTER': filter_val, 'samples': sample_gts,
            })
    return variants, samples


def is_fixed(gt_dict, sample_names):
    """Return True if all samples are homozygous alt (1/1)."""
    return len(sample_names) > 0 and all(
        gt_dict.get(s, './.') in ['1/1', '1|1'] for s in sample_names
    )


def main():
    parser = argparse.ArgumentParser(
        description='Variant annotation and classification (custom GFF+FASTA version)'
    )
    parser.add_argument('--vcf',    required=True, help='Input VCF file')
    parser.add_argument('--gff',    required=True, help='GFF3 annotation file')
    parser.add_argument('--ref',    required=True, help='Reference genome FASTA file')
    parser.add_argument('--outdir', default='./annotation_result', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("[1/4] Loading reference genome...")
    sequences = load_fasta_dict(args.ref)
    print(f"  Loaded {len(sequences)} sequences")

    print("[2/4] Parsing GFF file...")
    cds_intervals = parse_gff(args.gff)
    total_cds = sum(len(v) for v in cds_intervals.values())
    print(f"  Loaded {total_cds} CDS regions")

    print("[3/4] Parsing VCF file...")
    variants, samples = parse_vcf(args.vcf)
    print(f"  Loaded {len(variants)} variant sites, samples: {samples}")

    print("[4/4] Annotating variant effects...")
    mutations = {
        'missense':         [], 'missense_fixed':   [],
        'frameshift':       [], 'frameshift_fixed': [],
        'nonsense':         [], 'nonsense_fixed':   [],
    }

    in_cds_count = 0
    for i, var in enumerate(variants):
        cds_list = find_all_cds(var['CHROM'], var['POS'], cds_intervals)
        if not cds_list:
            continue

        in_cds_count += 1
        seen_effects = set()

        for cds_info in cds_list:
            effect, gene_id, mrna_id, aa_pos, (ref_aa, mut_aa) = predict_variant_effect(
                var['CHROM'], var['POS'], var['REF'], var['ALT'],
                cds_info, sequences
            )
            if not effect or effect in ('synonymous', 'unknown'):
                continue

            dedup_key = (var['CHROM'], var['POS'], effect, gene_id)
            if dedup_key in seen_effects:
                continue
            seen_effects.add(dedup_key)

            fixed = is_fixed(var['samples'], samples)
            entry = {
                'CHROM':       var['CHROM'],  'POS':        var['POS'],
                'REF':         var['REF'],    'ALT':        var['ALT'],
                'gene_id':     gene_id,       'mrna_id':    mrna_id,
                'effect':      effect,        'aa_position': aa_pos or '',
                'ref_aa':      ref_aa or '',  'mut_aa':     mut_aa or '',
                'is_fixed':    fixed,
            }
            for s in samples:
                entry[f'{s}_GT'] = var['samples'].get(s, './.')

            if effect == 'missense':
                mutations['missense'].append(entry)
                if fixed: mutations['missense_fixed'].append(entry)
            elif effect == 'frameshift':
                mutations['frameshift'].append(entry)
                if fixed: mutations['frameshift_fixed'].append(entry)
            elif effect == 'nonsense':
                mutations['nonsense'].append(entry)
                if fixed: mutations['nonsense_fixed'].append(entry)

        if (i + 1) % 5000 == 0:
            print(f"  Processed {i+1}/{len(variants)} variants, {in_cds_count} in CDS")

    print("\n=== Mutation Summary ===")
    for k, v in mutations.items():
        print(f"  {k}: {len(v)}")

    fieldnames_base = [
        'CHROM', 'POS', 'REF', 'ALT', 'gene_id', 'mrna_id',
        'effect', 'aa_position', 'ref_aa', 'mut_aa', 'is_fixed'
    ]
    gt_fields = [f'{s}_GT' for s in samples]
    fieldnames = fieldnames_base + gt_fields

    for category, var_list in mutations.items():
        output_file = os.path.join(args.outdir, f'{category}_mutations.csv')
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            if var_list:
                writer.writerows(var_list)
        print(f"  Written: {output_file} ({len(var_list)} sites)")

    print(f"\n[Done] Results saved to: {args.outdir}/")


if __name__ == '__main__':
    main()
