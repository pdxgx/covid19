#!/usr/bin/env python
"""
peptide_hits_plots.py

Generates tables of kmer coverage at each amino acid of a full viral peptide
sequence.

We ran

python peptide_hits_plots.py
    --fasta /path/to/viral/fasta --kmers /path/to/kmer/csv

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
    with open(args.kmers) as kmer_stream:
        kmer_stream.readline() # kill header
        for line in kmer_stream:
            allele, kmer = line.split(',')[:2]
            allele = allele[:5] # just get HLA-A, HLA-B, or HLA-C
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

