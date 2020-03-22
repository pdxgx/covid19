#!/usr/bin/env python3
"""
peptide_hits_plots.py

Generates tables of kmer coverage at each amino acid of a full viral peptide
sequence.

We ran

python3 peptide_hits_plots.py -f covid.fa -k covid_recode_conserved.csv -o covid
python3 peptide_hits_plots.py -f sars.fa -k covid_recode_conserved.csv -o sars

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
        '--conserved-kmers', '-k', type=str, required=True,
        help='kmer CSV')
    parser.add_argument(
        '--all-kmers', '-a', type=str, required=True,
        help='kmer CSV')
    parser.add_argument(
        '--out', '-o', type=str, required=True,
        help='output file basename'
        )
    args = parser.parse_args()

    id_to_name = {
            'YP_009725296.1' : 'ORF7b',
            'YP_009725295.1' : 'Last 12 AA of ORF1a',
            'YP_009725255.1' : 'ORF10',
            'YP_009724397.2' : 'N',
            'YP_009724396.1' : 'ORF8',
            'YP_009724395.1' : 'ORF7a',
            'YP_009724394.1' : 'ORF6',
            'YP_009724393.1' : 'M',
            'YP_009724392.1' : 'E',
            'YP_009724391.1' : 'ORF3a',
            'YP_009724390.1' : 'S',
            'YP_009724389.1' : 'ORF1ab'
        }
    name_order = ['Last 12 AA of ORF1a', 'ORF1ab', 'S', 'ORF3a', 'E', 'M',
                  'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF10']

    from collections import defaultdict
    kmer_sets = [defaultdict(set), defaultdict(set)]
    import csv
    for q, kmer_file in enumerate([args.conserved_kmers, args.all_kmers]):
        with open(kmer_file) as kmer_stream:
            kmer_reader = csv.reader(kmer_stream)
            for i, row in enumerate(kmer_reader):
                if not i:
                    if ('peptide' in row[0].lower()
                        and 'allele' in row[1].lower()):
                        reverse_order = True
                    elif ('peptide' in row[1].lower()
                          and 'allele' in row[0].lower()):
                        reverse_order = False
                    else:
                        raise RuntimeError(
                                'Input CSV does not have peptide/allele in '
                                'first and second columns'
                            )
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
                kmer_sets[q][allele].add(kmer)
                kmer_sets[q]['All'].add(kmer)

    import matplotlib.pyplot as plt
    #import seaborn as sns
    #sns.set()
    from pandas import DataFrame
    cov_dists_conserved = [
                defaultdict(lambda: [0 for i in range(len(viral_seq))])
                for p in range(len(name_order))
            ]
    cov_dists_all = [
                defaultdict(lambda: [0 for i in range(len(viral_seq))])
                for p in range(len(name_order))
            ]
    colors = { 'All' : '#787a64',
               'HLA-A' : '#1b9e77',
               'HLA-B' : '#d95f02',
               'HLA-C' : '#7570b3' }
    alpha_colors = { 'All' : '#ddded8',
                     'HLA-A' : '#c6e7dd',
                     'HLA-B' : '#f6d7c0',
                     'HLA-C' : '#dddbec' }
    with open(args.fasta) as fasta_stream:
        line = fasta_stream.readline().strip()
        if not line:
            raise RuntimeError('First line of FASTA is blank.')
        while True:
            assert line.startswith('>')
            contig = line[1:]
            for my_id in id_to_name:
                if my_id in contig:
                    p = name_order.index(id_to_name[my_id])
                    break
            viral_seq = []
            while True:
                line = fasta_stream.readline().strip()
                if line.startswith('>') or not line:
                    break
                viral_seq.append(line.strip())
            viral_seq = ''.join(viral_seq)
            for allele in ['All', 'HLA-A', 'HLA-B', 'HLA-C']:
                cov_dists_conserved[p][allele][0] = 0
                cov_dists_all[p][allele][0] = 0
            for kmer_size in [8, 9, 10, 11, 12]: # only these kmer sizes work
                for i in range(len(viral_seq) - kmer_size + 1):
                    for allele in colors:
                        if viral_seq[i:i+kmer_size] in kmer_sets[0][allele]:
                            for j in range(kmer_size):
                                cov_dists_conserved[p][allele][i+j] += 1
                        if viral_seq[i:i+kmer_size] in kmer_sets[1][allele]:
                            for j in range(kmer_size):
                                cov_dists_all[p][allele][i+j] += 1
            if not line:
                break
            continue
    plt.clf()
    all_cov_dists = defaultdict(list)
    for allele in colors:
        for p in range(len(cov_dists_all)):
            all_cov_dists[allele] += cov_dists_all[p][allele]
    for allele in all_cov_dists:
        with open(args.out + '.' + allele + '.full.csv', 'w') as all_stream:
            for el in all_cov_dists[allele]:
                print(el, file=all_stream)
    fig, axes = plt.subplots(
                        nrows=4, ncols=1,
                        sharex=True, sharey=True
                    )
    plt.setp(axes, xticks=range(0, len(all_cov_dists['All']), 2000),
                yticks=range(0, 50, 10))
    max_x = len(all_cov_dists['All']) + 10
    min_x = -10
    max_y = max(all_cov_dists['All']) + 2
    min_y = -0.2
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)
    plt.xlabel('Position')
    for k, allele in enumerate(colors):
        axes[k].vlines(x=range(len(all_cov_dists['All'])),
                       ymin=0, ymax=all_cov_dists[allele],
                       color=alpha_colors[allele])
        axes[k].text(max_x * 0.9, max_y * 1.03, allele)
    all_cov_dists = defaultdict(list)
    for allele in colors:
        for p in range(len(cov_dists_all)):
            all_cov_dists[allele] += cov_dists_conserved[p][allele]
    for allele in all_cov_dists:
        with open(args.out + '.' + allele + '.conserved.csv', 'w') as all_stream:
            for el in all_cov_dists[allele]:
                print(el, file=all_stream)
    for k, allele in enumerate(colors):
        axes[k].vlines(x=range(len(all_cov_dists['All'])),
                       ymin=0, ymax=all_cov_dists[allele],
                       color=colors[allele])
    fig.text(0.015, 0.5, 'Kmer count', ha='center', va='center',
             rotation='vertical')
    plt.tight_layout()

    plt.savefig(
            args.out + '.pdf'
        )
