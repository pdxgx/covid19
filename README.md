# covid19

[Read our paper in Journal of Virology](https://dx.doi.org/10.1128/JVI.00510-20)


## Software requirements

*R and R packages*

We ran `R` scripts with [`R` v3.6.2](https://ftp.osuosl.org/pub/cran/src/base/R-3/R-3.6.2.tar.gz). We used the following packages installable from CRAN:

* `dplyr` (v0.85)
* `rvest` (v0.3.5)
* `rworldmap` (v1.3-6)
* `ggplot2` (v3.3.0)

*Python and python packages:*

We ran python scripts with Python 3.7, installed using [Anaconda](https://www.anaconda.com/distribution/). See docstrings of `*.py` files for further information on reproducing results. Packages were installed with `pip3` and include:

* `matplotlib`
* `mmh3`

## Data

*SARS-CoV2 and other coronavirus protein sequences*

We manually downloaded the full polyprotein 1ab (ORF1ab), spike (S) protein, membrane (M) protein, envelope (E) protein, and nucleocapsid (N) protein sequence FASTA files for the following 34 coronavirus species from the [National Center of Biotechnology Information (NCBI) Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/). Note that for [SARS-CoV-2](https://www.ncbi.nlm.nih.gov/genome/?term=txid2697049) and [SARS-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid694009) we downloaded FASTA files for all additional annotated protein sequences comprising their full viral proteomes:


| Name | ORF1ab | Spike | Envelope | Membrane | Nucleocapsid |
| ---- | ------ | ----- | -------- | -------- | ------------ |
| [SARS-CoV-2\*](https://www.ncbi.nlm.nih.gov/genome/?term=txid2697049) |	YP_009724389.1 | YP_009724390.1 | YP_009724392.1 | YP_009724393.1 | YP_009724397.2 |
| [SARS-CoV\*](https://www.ncbi.nlm.nih.gov/genome/?term=txid694009) | NP_828849.2 | NP_828851.1 | NP_828854.1 | NP_828855.1 | NP_828858.1 |
| [OC43](https://www.ncbi.nlm.nih.gov/genome/?term=txid31631) | YP_009555238.1 | YP_009555241.1 | YP_009555243.1 | YP_009555244.1 | YP_009555245.1 |
| [Bovine-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid11128) | NP_150073.2 | NP_150077.1 | NP_150081.1 | NP_150082.1 | NP_150083.1 |
| [HKU24](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501960) | YP_009113022.1 | YP_009113025.1 | YP_009113028.1 | YP_009113029.1 | YP_009113031.1 |
| [HKU1](https://www.ncbi.nlm.nih.gov/genome/?term=txid290028) | YP_173236.1 | YP_173238.1 | YP_173240.1 | YP_173241.1 | YP_173242.1 |
| [MHV](https://www.ncbi.nlm.nih.gov/genome/?term=txid11138) | AAU06353.1 | AAU06356.1 | AAU06359.1 | AAU06360.1 | NP_045302.1 |
| [Rat-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid31632) | YP_003029844.1 | YP_003029848.1 | YP_003029850.1 | YP_003029851.1 | YP_003029852.1 |
| [Bat-BCoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501961) | YP_009072438.1 | YP_009072440.1 | YP_009072442.1 | YP_009072443.1 | YP_009072446.1 |
| [Hedgehog-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid1965093) | YP_009513008.1 | YP_009513010.1 | YP_009513016.1 | YP_009513017.1 | YP_009513018.1 |
| [MERS-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid1335626) | YP_009047202.1 | YP_009047204.1 | YP_009047209.1 | YP_009047210.1 | YP_009047211.1 |
| [HKU4](https://www.ncbi.nlm.nih.gov/genome/?term=txid694007) | YP_001039952.1 | YP_001039953.1 | YP_001039958.1 | YP_001039959.1 | YP_001039960.1 |
| [HKU5](https://www.ncbi.nlm.nih.gov/genome/?term=txid694008) | YP_001039961.1 | YP_001039962.1 | YP_001039967.1 | YP_001039968.1 | YP_001039969.1 |
| [GCCDC1](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501962) | YP_009273004.1 | YP_009273005.1 | YP_009273007.1 | YP_009273008.1 | YP_009273009.1 |
| [HKU9](https://www.ncbi.nlm.nih.gov/genome/?term=txid694006) | YP_001039970.1 | YP_001039971.1 | YP_001039973.1 | YP_001039974.1 | YP_001039975.1 |
| [HKU14](https://www.ncbi.nlm.nih.gov/genome/?term=txid1160968) | YP_005454239.1 | YP_005454245.1 | YP_005454247.1 | YP_005454248.1 | YP_005454249.1 |
| [CDPHE15](https://www.ncbi.nlm.nih.gov/genome/?term=txid1913643) | YP_008439200.1 | YP_008439202.1 | YP_008439204.1 | YP_008439205.1 | YP_008439206.1 |
| [HKU10](https://www.ncbi.nlm.nih.gov/genome/?term=txid1244203) | YP_006908641.2 | YP_006908642.1 | YP_006908644.1 | YP_006908645.1 | YP_006908646.1 |
| [BtRf-AlphaCoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501926) | YP_009199789.1 | YP_009199790.1 | YP_009199792.1 | YP_009199793.1 | YP_009199794.1 |
| [229E](https://www.ncbi.nlm.nih.gov/genome/?term=txid11137) | ARU07599.1 | ARU07601.1 | ARU07603.1 | ARU07604.1 | ARU07605.1 |
| [LuchengRn-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid1508224) | YP_009336483.1 | YP_009336484.1 | YP_009336485.1 | YP_009336486.1 | YP_009336487.1 |
| [Ferret-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid1264898) | YP_009256195.1 | YP_009256197.1 | YP_009256199.1 | YP_009256200.1 | YP_009256201.1 |
| [Mink-CoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid1913642) | YP_009019180.1 | YP_009019182.1 | YP_009019184.1 | YP_009019185.1 | YP_009019186.1 |
| [Bat-CoV-1A](https://www.ncbi.nlm.nih.gov/genome/?term=txid694000) | YP_001718603.1 | YP_001718605.1 | YP_001718607.1 | YP_001718608.1 | YP_001718609.1 |
| [HKU8](https://www.ncbi.nlm.nih.gov/genome/?term=txid694001) | YP_001718610.1 | YP_001718612.1 | YP_001718614.1 | YP_001718615.1 | YP_001718616.1 |
| [BtMr-AlphaCoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501927) | YP_009199608.1 | YP_009199609.1 | YP_009199611.1 | YP_009199612.1 | YP_009199613.1 |
| [BtNv-AlphaCoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501928) | YP_009201729.1 | YP_009201730.1 | YP_009201732.1 | YP_009201733.1 | YP_009201734.1 |
| [Porcine-EDV](https://www.ncbi.nlm.nih.gov/genome/?term=txid28295) | NP_598309.2 | NP_598310.1 | NP_598312.1 | NP_598313.1 | NP_598314.1 |
| [BtCoV512](https://www.ncbi.nlm.nih.gov/genome/?term=txid693999) | YP_001351683.1 | YP_001351684.1 | YP_001351686.1 | YP_001351687.1 | YP_001351688.1 |
| [HKU2](https://www.ncbi.nlm.nih.gov/genome/?term=txid693998) | YP_001552234.1 | YP_001552236.1 | YP_001552238.1 | YP_001552239.1 | YP_001552240.1 |
| [NL63](https://www.ncbi.nlm.nih.gov/genome/?term=txid277944) | YP_003766.2 | YP_003767.1 | YP_003769.1 | YP_003770.1 | YP_003771.1 |
| [NL63-related](https://www.ncbi.nlm.nih.gov/genome/?term=txid2501929) | YP_009328933.1 | YP_009328935.1 | YP_009328937.1 | YP_009328938.1 | YP_009328939.1 |
| [FCoV](https://www.ncbi.nlm.nih.gov/genome/?term=txid12663) | YP_004070193.2 | YP_004070194.1 | YP_004070197.1 | YP_004070198.1 | YP_004070199.1 |
| [TGE](https://www.ncbi.nlm.nih.gov/genome/?term=txid11149) | NP_058422.1 | NP_058424.1 | NP_058426.1 | NP_058427.2 | NP_058428.1 |

*HLA allele frequency data*

We designed and executed a custom `R` script (`HLA_frequencies.R`) to scrape all population and human leukocyte antigen (HLA) allele frequency from the [Allele Frequency Net Database](http://www.allelefrequencies.net/).


## Data analysis

*Protein sequence alignments*

For each protein class (i.e. ORF1ab, S, M, E, N), all 34 coronavirus sequences were aligned using the web-based [`Clustal Omega`](https://www.ebi.ac.uk/Tools/msa/clustalo/) multisequence aligner tool employing default parameters: sequence type [Protein], output alignment format [clustal_num], dealign [false], mBed-like clustering guide-tree [true], mBed-like clustering iteration [true], number of combined iterations [0], maximum guide tree iterations [-1], and maximum HMM iterations [-1].

*Binding affinity*

Using the above FASTA sequences, we used netchop v3.0 "C-term" model with a cleavage threshold of 0.1 to filter peptides that were not predicted to undergo canonical MHC class I antigen processing via proteasomal cleavage.

To assess binding affinity, we kmerized the protein FASTA into peptides of length 8-12.
``` sed 's/>.*/|/g' protein_sequence.fasta > protein_sequence_with_pipe.fasta ```
``` perl -pe 's/\s+//g' protein_sequence_with_pipe.fasta > protein_no_space.fasta ```

From the netchop v3.0 output, we further processed the file to get rid of NA lines and added row numbers as a column.
``` sed '/NA/d' netchop_out.txt > netchop_out_parsed.txt ``` 
``` awk -F'\t' 'NR>1{$0=$0"\t"NR-1} 1' netchop_out_parsed.txt > netchop_numbered.txt ```




