# deago
R package for generating HTML reports from differential expression analyis and GO term enrichment.

[![Build Status](https://travis-ci.org/sanger-pathogens/deago.svg?branch=master)](https://travis-ci.org/sanger-pathogens/deago)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/deago/blob/master/LICENSE)   
[![codecov](https://codecov.io/gh/sanger-pathogens/deago/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/deago)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [R](#r)
    * [Running the tests](#running-the-tests)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction
Generation of user-friendly HTML reports from differential expression analyis (DESeq2) and GO term enrichment (topGO) of count data using knitr. User provides count data and 'deago' will run the analysis and generate a HTML summary report containing QC plots, DE genes and top 30 GO terms.

## Installation
deago has the following dependencies:

### Required dependencies
* R >= 3.2
* devtools
* DESeq2
* topGO

If you encounter an issue when installing deago please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/deago/issues) or email us at path-help@sanger.ac.uk.

### R
Install the latest version of this package by entering the following in R:
```
install.packages("devtools")
library(devtools)
install_github("sanger-pathogens/deago")
```

### Running the tests
The test can be run from the top level directory, replacing <version> with the downloaded version:  
```
R CMD build --no-build-vignettes .
R CMD check deago_<version>.tar.gz --no-vignettes --as-cran --no-manual
```
## License
deago is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/deago/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/deago/issues) or email path-help@sanger.ac.uk.

## Citation
If you use this software please cite:

__DESeq2__:
__Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2__  
Love, M.I., Huber, W., Anders, S., Genome Biology 2014, 15:550. 10.1186/s13059-014-0550-8

__topGO__:
__topGO: Enrichment Analysis for Gene Ontology__  
Alexa, A. and Rahnenfuhrer, J., 2016, r BiocStyle::pkg_ver('topGO')
