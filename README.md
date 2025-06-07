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
- [Data Formatting & Required Files](#data-formatting--required-files)
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

### Overview

The following figure extracted shows a general overview of the overall EX-TERR pipeline:

### Data Formatting & Required Files

EX-TERR requires specific file formats and data preparation in order to execute pipeline.
Details and examples can be found after the summary table. 

| Required File | Brief Description |
|-----------------|-----------------|
| Genome-wide association studies (`.txt`)   | Summary statistics for external data from genome-wide association studies and  outcomes. "Outcome" summary statistics show the circular association between outcome variants and outcome phenotypes, as generated through REGENIE.   |
| Genotype Data <br> (`.bim/.bed/.fam`)  | Genotype data stored as PLINK binary files for outcome data.  Also includes files for allele freqency (.frq) and referenece allele. |
| No. of Blocks (`.txt`)   | A simple one-column file indicating the number of "sets" each genotype file is divided into, if formatted as such. |
| Traits lists (`.txt`) | List of traits corresponding to number of external GWAS and outcome traits tested. This pipeline can test multiple outcomes at once. |
| Corrections (`.txt`) | Cofactors of genotype information required for regression analyses, corresponding to participants' age, sex and desired number of genetic principal components (PCs). | 
| Phenotype (`.txt`) | Phenotypic information for outcome and validation steps. Also present in various formats | 
| Cross validation groups (`.txt`) | Pre-generated file with columns `<chr> <set> <ids>` corresponding to chromosome, set, and group identity (1-5) for the purpose of cross-validation. Default setting is for 5 folds. |
| Masking data information (`.txt`) | Pre-generated file signifying which GWAS traits to mask for each outcome. |
| Rotation information (`.txt`) | Pre-generated file which generate the rotation matrix which will be applied to genotype data, derived from training set genotype data. |


#### 1. Genome-wide association study summary statistics 
Genome-wide association study (GWAS) summary statistics represent the association between variants and a specific
trait or disease. The information from summary statistics used for the EX-TERR pipeline should not overlap with 
outcome phenotypes in order to avoid circularity and ensure unbiased effect estimates. 

The columns of the GWAS summary statistics follow the REGENIE output format, with the header:

<div align="center"> 
  
  `<rsid> <chr> <pos> <a0> <a1> <beta> <beta_se> <N>`
  
</div>

The header represents the following information:

<div align="center">
  
| Column | Description | 
|--------|-----------------|
| `rsid` | SNP identifier (Reference SNP ID) e.g. rs123456 |
| `chr` | Chromosome number (1-22) |
| `pos` | Base pair coordinate representing genomic position of the SNP |
| `a0` | Reference allele |
| `a1` | Alternate (effect) allele, which `beta` refers to |  
| `beta` | Regression coefficent, key input of EX-TERR. Estimated effect size `a1` on the trait/phenotype |
| `beta_se` | Standard error of the beta estimate | 
| `N` | Sample size used for the association at that variant | 

</div>

Note that it is important to identify which allele is the reference versus alternate allele for the study.
This can be further specified in the `ref_allele` file (further details 
found in the [Genotype data section](#2-genotype-data)).

#### 2. Genotype data

The genotypic data correponds to the group of individuals of which the outcome phenotypes are obtained from,
and used for validation within the EX-TERR pipeline. This data should be stored as PLINK binary files (`.bim`/`.fam`/`.bed`/) files as various sets. The following sets 
of participant should be available:

1. **Full cohort**
2. **Training set:** Used to train the MARS algorithm (e.g. 80\% of participants)
3. **Test set:** Used to verify MARS output (e.g. 20\% of participants)

Each set of genotypic data should be available cumulatively (all variants; `genotype_all`), as well as per chromosome 
(e.g. `genotype_chr1`, `genotype_chr2`...`genotype_chr22`). Chromosome data can be further divided into sets for 
larger chromosomes in neccesary (e.g. `genotype_chr1_set1` and `genotype_chr1_set2`). 
Note that a key aspect of EX-TERR is that it **does not** required an initial participant train-test split,
due to a downstream cross-validation step on the genotypic information. Thus, the pipeline can be trained on 
all available participants if desired. 

In addition to the core genotype information, details on allele information should also be prepared. 
Allele frequency information should be available for outcome variants file in the PLINK
format (`.frq`) , with the following column information:


Additionally, it is good practice to include `ref_allele` files to signify which allele from the genotype is
the alternate/effect allele. 



### Part 1: LDPRED2 Conversion

## License

## Contact Information
For questions, feedback or other inquiries regarding the EX-TERR pipeline, please contact 
Ann Le (annl.37@hotmail.com) or Guillaume Paré (pareg@mcmaster.ca). 

## Citation
