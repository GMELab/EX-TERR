# EX-TERR

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
- [Preparation](#preparation-data-formatting--getting-rotations)
- [Part 1: LDpred2 Conversion](#part-1-ldpred2-conversion)

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
This package requires multiple dependencies: `data.table`, `earth`, `caret`, `bigsnpr`, `magrittr`, `dplyr`, `pscl`, `fmsb`
and `glmnet`.

## Workflow

### Preparation: Data Formatting & Getting Rotations

EX-TERR requires specific file formats and data preparation in order to execute pipeline.
Summary of the files required can be seen below:

| Required File | Brief Description |
|-----------------|-----------------|
| Genome-wide association studies <br> (.txt)   | Summary statistics for external data from genome-wide association studies and UKB outcomes. UKB summary statistics show the association between UKB variants and UKB outcomes, as generated through REGENIE.   |
| Genotype Data <br> (bim/bed/fam)  | Genotype data stored as PLINK binary files for outcome data.  Also includes files for allele freqency (.frq) and referenece allele. |
| No. of Blocks   | A simple one-column file indicating the number of "sets" each genotype file is divided into, if formatted as such. |
| Traits lists <br> (.txt) | List of traits corresponding to number of external GWAS and outcome traits tested. This pipeline can test multiple outcomes at once. |
| Corrections <br> (.txt) | Cofactors of genotype information required for regression analyses, corresponding to participants' age, sex and desired number of genetic principal components (PCs). | 
| Phenotype <br> (.txt) | Phenotypic information for outcome and validation steps. Also present in various formats | 
| Cross validation groups | Pre-generated file with columns `<chr> <set> <ids>` corresponding to chromosome, set, and group identity (1-5) for the purpose of cross-validation. Default setting is for 5 folds. |

### Part 1: LDPRED2 Conversion

## Contact Information
For questions, feedback or other inquiries regarding the EX-TERR pipeline, please contact 
Ann Le (annl.37@hotmail.com) or Guillaume Paré (pareg@mcmaster.ca). 
