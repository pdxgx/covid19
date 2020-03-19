#!/usr/bin/env python3
"""
peptide_hits_plots.py

Generates tables of kmer coverage at each amino acid of a full viral peptide
sequence.

We ran

python3 peptide_hits_plots.py -f YP_009724389.1.fa -k hla_cleavedflurry_filt_50.csv -o covid19.flurry
python3 peptide_hits_plots.py -f YP_009724389.1.fa -k hla_cleavednugs_filt_50.csv -o covid19.nugs
python3 peptide_hits_plots.py -f NP_828849.2.fa -k sars_hla_cleavedflurry_filt_50.csv -o sars.flurry
python3 peptide_hits_plots.py -f NP_828849.2.fa -k sars_hla_cleavednugs_filt_50.csv -o sars.nugs

kmer CSV has columns allele, peptide, ic50, less_than_50, occ, but only the
columns allele and peptide are used.

Outputs coverage at each position of input peptide for each of HLA-A, HLA-B,
HLA-C, and all.
"""

import seaborn

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

    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set()
    from pandas import DataFrame
    with open(args.fasta) as fasta_stream:
        line = fasta_stream.readline().strip()
        if not line:
            raise RuntimeError('First line of FASTA is blank.')
        while True:
            assert line.startswith('>')
            contig = line[1:]
            viral_seq = []
            while True:
                line = fasta_stream.readline().strip()
                if line.startswith('>') or not line:
                    break
                viral_seq.append(line.strip())
            viral_seq = ''.join(viral_seq)
            cov_dists = defaultdict(lambda: [0 for i in range(len(viral_seq))])
            for kmer_size in [8, 9, 10, 11, 12]: # only these kmer sizes work
                for i in range(len(viral_seq) - kmer_size + 1):
                    for allele in kmer_sets:
                        if viral_seq[i:i+kmer_size] in kmer_sets[allele]:
                            for j in range(kmer_size):
                                cov_dists[allele][i+j] += 1

            plt.clf()
            fig, axes = plt.subplots(
                                nrows=len(cov_dists), ncols=1,
                                sharex=True, sharey=True
                            )
            plt.setp(axes, xticks=range(0, len(viral_seq), 1000),
                        yticks=range(0, 50, 5))
            max_x = len(viral_seq) + 10
            min_x = -10
            max_y = max(cov_dists['all']) + 1
            min_y = -0.2
            plt.xlim(min_x, max_x)
            plt.ylim(min_y, max_y)
            colors = { 'all' : '#787a64',
                       'HLA-A' : '#1b9e77',
                       'HLA-B' : '#d95f02',
                       'HLA-C' : '#7570B3' }
            for k, allele in enumerate(colors.keys()):
                with open(args.out + '.' + contig + '.' + allele + '.csv',
                          'w') as cov_stream:
                    for i in range(len(cov_dists[allele])):
                        print(cov_dists[allele][i], file=cov_stream)
                axes[k].vlines(x=range(len(viral_seq)),
                           ymin=0, ymax=cov_dists[allele],
                           color=colors[allele])
                axes[k].text(max_x * 0.9, max_y * 1.03, allele)
            plt.xlabel('position')
            fig.text(0.015, 0.5, 'kmer count', ha='center', va='center',
                     rotation='vertical')
            plt.tight_layout()
            plt.savefig(
                    args.out + '.' + contig + '.pdf'
                )
            if not line:
                break
            continue
