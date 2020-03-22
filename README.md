# covid19

## Software requirements

*R and R packages*

We ran `R` scripts with [`R` v3.6.2](https://ftp.osuosl.org/pub/cran/src/base/R-3/R-3.6.2.tar.gz). We used the following packages installable from CRAN:

* `dplyr` (v0.85)
* `rvest` (v0.3.5)
* `rworldmap` (v1.3-6)
* `ggplot2` (v3.3.0)

We also used [`R studio`(https://rstudio.com/products/rstudio/download/#download) for interactively generating figures and analyzing results.

## Data

*SARS-CoV2 and other coronavirus protein sequences*


*HLA frequency data*

Allele search in allelefrequencies.com

We scraped population and allele frequency data from allelefrequency.com. We used the `rvest` R package to parse the scraped html code, and the `rglobalmaps`.
