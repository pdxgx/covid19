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

    import seaborn as sns
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

            for allele in cov_dists:
                with open(args.out + '.' + contig + '.' + allele + '.csv',
                          'w') as cov_stream:
                    for i in range(len(cov_dists[allele])):
                        print(cov_dists[allele][i], file=cov_stream)
                        sns.barplot(x = 'position', y = 'kmer count',
                                    data = cov_dists[allele],
                                    palette = 'hls')
                        plt.savefig()
            if not line:
                break
            continue
