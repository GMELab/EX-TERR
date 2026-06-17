#!/bin/sh

#' @export
baseline_prs_asso <- function(mode, age, sex, pcs, input_dir, output_dir = NULL) {
        suppressMessages(library("data.table"))

        if (!is.null(output_dir)) {
                dir.create(output_dir, recursive = TRUE)
        }
        if (!mode %in% c("default", "no_mask", "gwas", "with_name")) {
                stop("mode should be either 'default' or 'no_mask' or 'gwas' or 'with_name'")
        }

        standardization <- function(x) {
                return((x - mean(x)) / sd(x))
        }

        get_cont_asso <- function(trait, y, PRS) {
                pheno_resid <- resid(lm(y ~ age + sex + pcs))
                pheno_norm <- standardization(pheno_resid)
                adj_r2 <- 0
                for (j in seq_len(dim(PRS)[2]))
                {
                        temp <- summary(lm(pheno_norm ~ standardization(PRS[, j])))
                        if (temp$adj.r.squared > adj_r2) {
                                if (mode == "gwas" || mode == "with_name") {
                                        item <- c(trait, colnames(PRS)[j], coef(temp)[2, c(1, 2, 4)], temp$adj.r.squared)
                                } else {
                                        item <- c(trait, coef(temp)[2, c(1, 2, 4)], temp$adj.r.squared)
                                }
                                adj_r2 <- temp$adj.r.squared
                        }
                }
                return(item)
        }

        get_dicho_asso <- function(trait, y, PRS) {
                OR <- 0
                for (j in seq_len(dim(PRS)[2]))
                {
                        temp <- summary(glm(y ~ standardization(PRS[, j]) + age + sex + pcs, family = "binomial"))
                        if (exp(coef(temp)[2, 1]) > OR) {
                                if (mode == "gwas" || mode == "with_name") {
                                        item <- c(trait, colnames(PRS)[j], coef(temp)[2, c(1, 2, 4)], exp(coef(temp)[2, 1]))
                                } else {
                                        item <- c(trait, coef(temp)[2, c(1, 2, 4)], exp(coef(temp)[2, 1]))
                                }
                                OR <- exp(coef(temp)[2, 1])
                        }
                }
                return(item)
        }

        blocks_orig <- as.matrix(fread(file.path(input_dir, "Cross_validation_groups.txt"))) # Randomly divided into 5 groups
        if (mode == "gwas" || mode == "with_name") {
                testing_traits <- as.matrix(fread(file.path(input_dir, "mask_outcome_list_final.txt"), header = F))
        } else {
                testing_traits <- as.matrix(fread(file.path(input_dir, "outcome_list_final.txt"), header = F)) # 71 test traits
        }

        # Test association in discovery set
        Age <- as.matrix(fread(file.path(input_dir, "Age_disc.txt"), header = T))
        Sex <- as.matrix(fread(file.path(input_dir, "Sex_disc.txt"), header = T))
        PCs <- as.matrix(fread(file.path(input_dir, "PCs_disc.txt"), header = T))

        age <- Age[, 2]
        sex <- Sex[, 2]
        pcs <- PCs[, 2:11]

        # Get PRS_disc and PRS_val
        if (!file.exists(file.path(input_dir, "Baseline_PRS_disc.txt")) || !file.exists(file.path(input_dir, "Baseline_PRS_val.txt"))) {
                prs_disc <- as.matrix(fread(file.path(input_dir, "PRS_chr_1_1_1_disc.txt")))
                PRS_disc <- matrix(0, ncol = dim(prs_disc)[2], nrow = dim(prs_disc)[1])

                prs_val <- as.matrix(fread(file.path(input_dir, "PRS_chr_1_1_1_val.txt")))
                PRS_val <- matrix(0, ncol = dim(prs_val)[2], nrow = dim(prs_val)[1])

                for (i in seq_len(dim(blocks_orig)[1]))
                {
                        item <- blocks_orig[i, ]
                        prs_disc <- as.matrix(fread(file.path(input_dir, paste0("PRS_chr_", item[1], "_", item[2], "_", item[3], "_disc.txt"))))
                        PRS_disc <- PRS_disc + prs_disc

                        prs_val <- as.matrix(fread(file.path(input_dir, paste0("PRS_chr_", item[1], "_", item[2], "_", item[3], "_val.txt"))))
                        PRS_val <- PRS_val + prs_val
                }

                write.table(PRS_disc, file.path(input_dir, "Baseline_PRS_disc.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
                write.table(PRS_val, file.path(input_dir, "Baseline_PRS_val.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
        } else {
                PRS_disc <- as.matrix(fread(file.path(input_dir, "Baseline_PRS_disc.txt")))
                PRS_val <- as.matrix(fread(file.path(input_dir, "Baseline_PRS_val.txt")))
        }

        if (mode == "gwas") {
                Results_cont <- c("Traits", "Best_GWAS", "Beta", "beta_se", "pval", "adj_r2")
                Results_dicho <- c("Traits", "Best_GWAS", "Beta", "beta_se", "pval", "OR")
        } else if (mode == "with_name") {
                Results_cont <- c("Traits", "GWAS_name", "Beta", "beta_se", "pval", "adj_r2")
                Results_dicho <- c("Traits", "GWAS_name", "Beta", "beta_se", "pval", "OR")
        } else {
                Results_cont <- c("Traits", "Beta", "beta_se", "pval", "adj_r2")
                Results_dicho <- c("Traits", "Beta", "beta_se", "pval", "OR")
        }

        for (i in seq_len(length(testing_traits)))
        {
                trait <- testing_traits[i]
                phenos <- as.matrix(fread(file.path(input_dir, paste0(trait, "_disc.txt")), header = T))


                if (mode != "no_mask" && file.exists(file.path(input_dir, paste0(trait, "_masked.txt")))) {
                        mask_gwas <- as.matrix(fread(file.path(input_dir, paste0(trait, "_masked.txt"), header = F)))
                        mask_ids <- match(mask_gwas, colnames(PRS_disc))

                        PRS_disc_masked <- PRS_disc[, -mask_ids]
                } else {
                        PRS_disc_masked <- PRS_disc
                }

                y <- phenos[, 3, drop = F]
                if (max(y) == 1 && min(y) == 0) # dichotomous mode
                        {
                                item <- get_dicho_asso(trait, y, PRS_disc_masked)
                                Results_dicho <- rbind(Results_dicho, item)
                        } else # continuous mode
                {
                        item <- get_cont_asso(trait, y, PRS_disc_masked)
                        Results_cont <- rbind(Results_cont, item)
                }
        }

        if (!is.null(output_dir)) {
                write.table(Results_dicho, file.path(output_dir, "Baseline_PRS_regression_dichotomous_disc.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
                write.table(Results_cont, file.path(output_dir, "Baseline_PRS_regression_continuous_disc.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        }

        # Test association in validation set
        Age <- as.matrix(fread(file.path(input_dir, "Age_val.txt"), header = T))
        Sex <- as.matrix(fread(file.path(input_dir, "Sex_val.txt"), header = T))
        PCs <- as.matrix(fread(file.path(input_dir, "PCs_val.txt"), header = T))

        age <- Age[, 2]
        sex <- Sex[, 2]
        pcs <- PCs[, 2:11]

        if (mode == "gwas") {
                Results_cont <- c("Traits", "Best_GWAS", "Beta", "beta_se", "pval", "adj_r2")
                Results_dicho <- c("Traits", "Best_GWAS", "Beta", "beta_se", "pval", "OR")
        } else if (mode == "with_name") {
                Results_cont <- c("Traits", "GWAS_name", "Beta", "beta_se", "pval", "adj_r2")
                Results_dicho <- c("Traits", "GWAS_name", "Beta", "beta_se", "pval", "OR")
        } else {
                Results_cont_val <- c("Traits", "Beta", "beta_se", "pval", "adj_r2")
                Results_dicho_val <- c("Traits", "Beta", "beta_se", "pval", "OR")
        }

        for (i in seq_len(length(testing_traits)))
        {
                trait <- testing_traits[i]
                # phenos = as.matrix(fread(paste0("Traits_UKB/", trait, "/Pheno/", trait, "_val_updated.txt"), header = T))
                phenos <- as.matrix(fread(file.path(input_dir, paste0(trait, "_val.txt")), header = T))

                # if(file.exists(paste0("Masked_data/", trait, "_masked.txt")))
                if (file.exists(file.path(input_dir, paste0(trait, "_masked.txt"))) && mode != "no_mask") {
                        mask_gwas <- as.matrix(fread(file.path(input_dir, paste0(trait, "_masked.txt")), header = F))
                        mask_ids <- match(mask_gwas, colnames(PRS_disc))

                        PRS_val_masked <- PRS_val[, -mask_ids]
                } else {
                        PRS_val_masked <- PRS_val
                }

                y <- phenos[, 3, drop = F]
                if (max(y) == 1 && min(y) == 0) # dichotomous mode
                        {
                                item <- get_dicho_asso(trait, y, PRS_val_masked)
                                Results_dicho_val <- rbind(Results_dicho_val, item)
                        } else # continuous mode
                {
                        item <- get_cont_asso(trait, y, PRS_val_masked)
                        Results_cont_val <- rbind(Results_cont_val, item)
                }
        }

        if (!is.null(output_dir)) {
                write.table(Results_dicho_val, file.path(output_dir, "Baseline_PRS_regression_dichotomous_val.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
                write.table(Results_cont_val, file.path(output_dir, "Baseline_PRS_regression_continuous_val.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
        }

        return(list(
                dicho = Results_dicho, cont = Results_cont,
                dicho_val = Results_dicho_val, cont_val = Results_cont_val
        ))
}
