#!/bin/sh

#' @param mask Specifies whether or not to applying masking condition (leaving out most highly associated trait)
#' @param corrections_dir Path to directory containing corrections for Age, Sex and PCs. Files should be in the form <correction>_<flag>.txt
#' @param pheno_dir Path to directory containing phenotype information separated by trait and disc and val. Files are named <trait>_disc.txt and <trait>_val.txt
#' @return Returns results as element "results"
#' @export
multi_PRS_asso <- function(
    mask,
    corrections_dir,
    pheno_dir,
    prs_dir) {
        suppressMessages(library("data.table"))

        if (mask != "mask" && mask != "no_mask") {
                stop("Mask must be either 'mask' or 'no_mask'.")
        }
        if (mask == "no_mask") {
                flag_2 <- "_no_mask"
        } else {
                flag_2 <- ""
        }

        standardization <- function(x) {
                return((x - mean(x)) / sd(x))
        }

        Age <- as.matrix(fread(file.path(corrections_dir, "Age_disc.txt"), header = T))
        Sex <- as.matrix(fread(file.path(corrections_dir, "Sex_disc.txt"), header = T))
        PCs <- as.matrix(fread(file.path(corrections_dir, "PCs_disc.txt"), header = T))

        age <- Age[, 2]
        sex <- Sex[, 2]
        pcs <- PCs[, 2:11]

        # For association in discovery set
        one_dicho_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_dicho_disc_1", flag_2, ".txt")))) # This is one out group in discovery set
        one_cont_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_continuous_disc_1", flag_2, ".txt")))) # This is one out group in discovery set
        print(dim(one_dicho_prs))
        print(dim(one_cont_prs))

        one_dicho_PRS <- matrix(0, ncol = dim(one_dicho_prs)[2], nrow = dim(one_dicho_prs)[1])
        one_cont_PRS <- matrix(0, ncol = dim(one_cont_prs)[2], nrow = dim(one_cont_prs)[1])

        for (ids in 1:5)
        {
                one_cont_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_continuous_disc_", ids, flag_2, ".txt"))))
                one_cont_PRS <- one_cont_PRS + one_cont_prs

                one_dicho_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_dicho_disc_", ids, flag_2, ".txt"))))
                one_dicho_PRS <- one_dicho_PRS + one_dicho_prs
        }

        if (nrow(one_cont_PRS) > 0) {
                Results <- c("Traits", "Beta", "beta_se", "pval", "adj_r2")
                cont_trait <- colnames(one_cont_PRS)
                for (i in seq_along(cont_trait))
                {
                        trait <- cont_trait[i]

                        phenos <- as.matrix(fread(paste0(pheno_dir, "/", trait, "_disc.txt"), header = T))
                        pheno_resid <- resid(lm(phenos[, 3] ~ age + sex + pcs))
                        pheno_norm <- standardization(pheno_resid)

                        multi_PRS <- one_cont_PRS[, i]

                        temp <- summary(lm(pheno_norm ~ standardization(multi_PRS)))
                        item <- c(trait, coef(temp)[2, c(1, 2, 4)], temp$adj.r.squared)
                        Results <- rbind(Results, item)
                }
                write.table(Results, file.path(prs_dir, paste0("R7_LASSO_Asso_cont_disc", flag_2, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")
        }

        if (nrow(one_dicho_PRS) > 0) {
                Results <- c("Traits", "Beta", "beta_se", "pval", "OR")
                dicho_trait <- colnames(one_dicho_PRS)

                for (i in seq_along(dicho_trait))
                {
                        trait <- dicho_trait[i]
                        phenos <- as.matrix(fread(file.path(pheno_dir, paste0(trait, "_disc.txt")), header = T))
                        pheno_norm <- phenos[, 3, drop = F]
                        pheno_norm <- ifelse(pheno_norm != 0, 1, 0)

                        multi_PRS <- one_dicho_PRS[, i]

                        temp <- summary(glm(pheno_norm ~ standardization(multi_PRS) + age + sex + pcs, family = "binomial"))
                        item <- c(trait, coef(temp)[2, c(1, 2, 4)], exp(coef(temp)[2, 1]))
                        Results <- rbind(Results, item)
                }
                write.table(Results, file.path(prs_dir, paste0("R7_LASSO_Asso_dicho_disc", flag_2, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")
        }


        # For association in validation set
        one_dicho_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_dicho_val_1", flag_2, ".txt")))) # This is one out group in validation set
        one_cont_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_continuous_val_1", flag_2, ".txt")))) # This is one out group in validation set

        one_dicho_PRS <- matrix(0, ncol = dim(one_dicho_prs)[2], nrow = dim(one_dicho_prs)[1])
        one_cont_PRS <- matrix(0, ncol = dim(one_cont_prs)[2], nrow = dim(one_cont_prs)[1])

        for (ids in 1:5)
        {
                one_cont_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_continuous_val_", ids, flag_2, ".txt"))))
                one_cont_PRS <- one_cont_PRS + one_cont_prs

                one_dicho_prs <- as.matrix(fread(file.path(prs_dir, paste0("R6_LASSO_dicho_val_", ids, flag_2, ".txt"))))
                one_dicho_PRS <- one_dicho_PRS + one_dicho_prs
        }

        Age <- as.matrix(fread(file.path(corrections_dir, "Age_val.txt"), header = T))
        Sex <- as.matrix(fread(file.path(corrections_dir, "Sex_val.txt"), header = T))
        PCs <- as.matrix(fread(file.path(corrections_dir, "PCs_val.txt"), header = T))

        age <- Age[, 2]
        sex <- Sex[, 2]
        pcs <- PCs[, 2:11]

        if (nrow(one_cont_PRS) > 0) {
                Results <- c("Traits", "Beta", "beta_se", "pval", "adj_r2")
                cont_trait <- colnames(one_cont_PRS)
                for (i in seq_along(cont_trait))
                {
                        trait <- cont_trait[i]

                        phenos <- as.matrix(fread(file.path(pheno_dir, paste0(trait, "_val.txt")), header = T))
                        pheno_resid <- resid(lm(phenos[, 3] ~ age + sex + pcs))
                        pheno_norm <- standardization(pheno_resid)

                        multi_PRS <- one_cont_PRS[, i]

                        temp <- summary(lm(pheno_norm ~ standardization(multi_PRS)))
                        item <- c(trait, coef(temp)[2, c(1, 2, 4)], temp$adj.r.squared)
                        Results <- rbind(Results, item)
                }
                write.table(Results, file.path(prs_dir, paste0("R7_LASSO_Asso_cont_val", flag_2, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")
        }

        if (nrow(one_dicho_PRS) > 0) {
                Results <- c("Traits", "Beta", "beta_se", "pval", "OR")
                dicho_trait <- colnames(one_dicho_PRS)
                for (i in seq_along(dicho_trait))
                {
                        trait <- dicho_trait[i]

                        phenos <- as.matrix(fread(file.path(pheno_dir, paste0(trait, "_val.txt")), header = T))

                        pheno_norm <- phenos[, 3, drop = F]
                        pheno_norm <- ifelse(pheno_norm != 0, 1, 0)

                        multi_PRS <- one_dicho_PRS[, i]

                        temp <- summary(glm(pheno_norm ~ standardization(multi_PRS) + age + sex + pcs, family = "binomial"))
                        item <- c(trait, coef(temp)[2, c(1, 2, 4)], exp(coef(temp)[2, 1]))
                        Results <- rbind(Results, item)
                }
                write.table(Results, file.path(prs_dir, paste0("R7_LASSO_Asso_dicho_val", flag_2, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")
        }

        return(list(results = Results))
}
