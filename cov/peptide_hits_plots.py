#!/usr/bin/env python3
"""
peptide_hits_plots.py

Generates tables of kmer coverage at each amino acid of a full viral peptide
sequence.

We ran

pypy3 peptide_hits_plots.py -f YP_009724389.1.fa -k hla_cleavednugs_filt.csv -o covid19.nugs
pypy3 peptide_hits_plots.py -f YP_009724389.1.fa -k hla_cleavedflurry_filt.csv -o covid19.flurry

kmer CSV has columns allele, peptide, ic50, less_than_50, occ, but only the
columns allele and peptide are used.

Outputs coverage at each position of input peptide for each of HLA-A, HLA-B,
HLA-C, and all.
"""

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--fasta', '-f', type=str, required=True,
        help='fasta with viral peptide sequence')
    parser.add_argument(
        '--kmers', '-k', type=str, required=True,
        help='kmer CSV')
    parser.add_argument(
        '--out', '-o', type=str, required=True,
        help='output file basename'
        )
    args = parser.parse_args()


    viral_seq = []
    with open(args.fasta) as fasta_stream:
        fasta_stream.readline() # kill header
        for line in fasta_stream:
            viral_seq.append(line.strip())
    viral_seq = ''.join(viral_seq)

    from collections import defaultdict
    kmer_sets = defaultdict(set)
    import csv
    with open(args.kmers) as kmer_stream:
        kmer_reader = csv.reader(kmer_stream)
        for i, row in enumerate(kmer_reader):
            if not i:
                if 'peptide' in row[0] and 'allele' in row[1]:
                    reverse_order = True
                elif 'peptide' in row[1] and 'allele' in row[0]:
                    reverse_order = False
                else:
                    raise RuntimeError('Input CSV does not have peptide/allele in '
                                       'first and second columns')
            if not reverse_order:
                allele, kmer = row[:2]
            else:
                kmer, allele = row[:2]
            # just get HLA-A, HLA-B, or HLA-C
            if 'HLA-A' in allele:
                allele = 'HLA-A'
            elif 'HLA-B' in allele:
                allele = 'HLA-B'
            elif 'HLA-C' in allele:
                allele = 'HLA-C'
            kmer_sets[allele].add(kmer)
            kmer_sets['all'].add(kmer)

    cov_dists = defaultdict(lambda: [0 for i in range(len(viral_seq))])
    for kmer_size in [8, 9, 10, 11, 12]: # NOTE only these kmer sizes work
        for i in range(len(viral_seq) - kmer_size + 1):
            for allele in kmer_sets:
                if viral_seq[i:i+kmer_size] in kmer_sets[allele]:
                    for j in range(kmer_size):
                        cov_dists[allele][i+j] += 1

    for allele in cov_dists:
        with open(args.out + '.' + allele + '.csv', 'w') as cov_stream:
            for i in range(len(cov_dists[allele])):
                print(cov_dists[allele][i], file=cov_stream)
