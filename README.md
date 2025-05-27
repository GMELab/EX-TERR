 EX-TERR

> EX-TERR is a methodology for polygenic risk score (PRS) generation, utilizing supervised machine learning to 
combine and priortize polygenic genetic information from multiple different traits from related risk factors.
> Specifically, the pipeline utilizes multivariate adaptive regression splines (MARS) to recalibrate
> and determine the association strength of each input trait. EX-TERR applies MARS to multiple rotated
> genetic variant matrices generated through principal component analysis (PCA) techniques. Uniquely,
> EX-TERR conducts a 5-fold cross-validation (CV) at the genotype level, eliminating the need
> for an initial participant-level train-test split. This addresses potential power reductions
> should GWAS information be unavailable. Further details regarding development and further details
> of EX-TERR can be referred to in the [corresponding manuscript](https://example.com).
>
> This R package automates the pipeline using functions: 

## Table of Contents
- [Installation](#installation)
- [Contact Information](#contact-information)

## Installation
The direct installation of EX-TERR can be performed using `devtools::install_github`:

```sh
devtools::install_github("GMELab/EXTERR")
```

If issues are encountered with the above command, a local download can be performed
using the `git clone` command as follows:

```sh
# Enter target directory to download EX-TERR
cd directory

# Initialize git repository
git init

# Download locally using the git clone command
git clone https://github.com/GMELab/EX-TERR

```

EX-TERR can be loaded in `R` using `devtools::load_all`:

```R
setwd("/directory/R")
devtools::load_all(".")
```

## Workflow

## Contact Information
For questions, feedback or other inquiries regarding the EX-TERR pipeline, please contact 
Ann Le (annl.37@hotmail.com) or Guillaume Paré (pareg@mcmaster.ca). 
