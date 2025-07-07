#' @param flag Type of set: disc (discovery) or val (validation)
#' @param PC_std_threshold PC SD threshold for filtering of rotated genotype data
#' @param mask Specifies whether or not to applying masking condition (leaving out most highly associated trait)
#' @param corrections_dir Path to directory containing corrections for Age, Sex and PCs. Files should be in the form <correction>_<flag>.txt
#' @param pheno_dir Path to directory containing phenotype file. Name should be in the form Pheno_<flag>.txt
#' @param genotype_dir Path to directory containing genotype files. Directory contains directories Geno_<flag> and files <outcome_db>_final.fam
#' @return A list containing two elements: earth_cont & earth_dicho
#' @export
get_assoc <- function(flag,
                      PC_std_threshold,
                      mask = c("mask", "no_mask"),
                      corrections_dir,
                      genotype_dir,
                      pheno_dir,
                      outcome_db,
                      geno_disc_dir,
                      output_dir) {
  suppressMessages(library("data.table"))

  output_dir <- file.path(output_dir)
  if (mask == "mask") {
    flag_2 <- ""
  } else {
    flag_2 <- "_no_mask"
  }

  standardization <- function(x) {
    return((x - mean(x)) / sd(x))
  }

  standardizing_trait <- function(pheno_data) {
    pheno_norm <- as.numeric()
    for (i in seq_len(dim(pheno_data)[2]))
    {
      print("STD", i)
      if (max(pheno_data[, i], na.rm = T) == 1 && min(pheno_data[, i], na.rm = T) == 0) {
        Y_norm <- pheno_data[, i]
        Y_norm[which(is.na(Y_norm))] <- 0
      } else {
        pheno_data[which(is.na(pheno_data[, i])), i] <- mean(pheno_data[, i], na.rm = T)
        Y_norm <- standardization(pheno_data[, i])
      }
      pheno_norm <- cbind(pheno_norm, Y_norm)
    }
    colnames(pheno_norm) <- colnames(pheno_data)

    return(pheno_norm)
  }

  # Main function

  Age <- as.matrix(fread(file.path(corrections_dir, paste0("Age_", flag, ".txt")), header = T))
  Sex <- as.matrix(fread(file.path(corrections_dir, paste0("Sex_", flag, ".txt")), header = T))
  PCs <- as.matrix(fread(file.path(corrections_dir, paste0("PCs_", flag, ".txt")), header = T))

  phenos <- data.frame(fread(file.path(pheno_dir, paste0("Pheno_", flag, ".txt")), header = T))
  # pheno_data <- phenos[, -1, drop = FALSE]
  # pheno_norm <- standardizing_trait(pheno_data)

  # Match fam order
  fam <- data.frame(fread(file.path(genotype_dir, paste0("Geno_", flag), paste0(outcome_db, "_final.fam"))))
  phenos <- as.matrix(phenos[match(fam[, 1], phenos[, 1]), ])
  pheno_data <- phenos[, -1, drop = FALSE]
  pheno_norm <- standardizing_trait(pheno_data)

  earth_PRS <- list()
  for (ids in 1:5)
  {
    temp <- as.matrix(fread(file.path(geno_disc_dir, paste0("Earth_PRS_PC_std_threshold_", PC_std_threshold, "_index_", ids, "_", flag, flag_2, ".txt"))))
    temp[which(is.na(temp))] <- 0
    earth_PRS[[ids]] <- temp
  }

  earth_PRS[[6]] <- matrix(0, ncol = dim(temp)[2], nrow = dim(temp)[1])
  for (ids in 1:5)
  {
    earth_PRS[[6]] <- earth_PRS[[6]] + earth_PRS[[ids]]
  }

  Earth_NAs <- c("Trait", "cross_id")
  ids <- 6
  Earth_cont <- c("Trait", "beta", "beta_SE", "pval", "adj_r2")
  Earth_dicho <- c("Trait", "beta", "beta_SE", "pval", "OR")
  for (i in seq_len(dim(pheno_norm)[2]))
  {
    pheno <- colnames(pheno_norm)[i]
    print(pheno, i)
    if (max(pheno_data[, i]) == 1 && min(pheno_data[, i]) == 0) {
      if (max(earth_PRS[[ids]][, i]) != 0 && min(earth_PRS[[ids]][, i]) != 0) {
        temp3 <- summary(glm(pheno_norm[, i] ~ standardization(earth_PRS[[ids]][, i]) + Age[, 2] + Sex[, 2] + PCs[, c(2:11)], family = binomial))
        item3 <- c(pheno, coef(temp3)[2, c(1, 2, 4)], exp(coef(temp3)[2, 1]))
        Earth_dicho <- rbind(Earth_dicho, item3)
      } else {
        item3 <- c(pheno, "NA", "NA", "NA", "NA")
        Earth_dicho <- rbind(Earth_dicho, item3)
        item <- c(pheno, ids)
        Earth_NAs <- rbind(Earth_NAs, item)
      }
    } else {
      if (max(earth_PRS[[ids]][, i]) != 0 && min(earth_PRS[[ids]][, i]) != 0) {
        temp1 <- summary(lm(pheno_norm[, i] ~ standardization(earth_PRS[[ids]][, i]) + Age[, 2] + Sex[, 2] + PCs[, c(2:11)]))
        item1 <- c(pheno, coef(temp1)[2, c(1, 2, 4)], temp1$adj.r.squared)
        Earth_cont <- rbind(Earth_cont, item1)
      } else {
        item1 <- c(pheno, "NA", "NA", "NA", "NA")
        Earth_cont <- rbind(Earth_cont, item1)
        item <- c(pheno, ids)
        Earth_NAs <- rbind(Earth_NAs, item)
      }
    }
  }

  write.table(Earth_cont, file.path(output_dir, paste0("Earth_cont_PC_std_threshold_", PC_std_threshold, "_", ids, "_", flag, flag_2, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(Earth_dicho, file.path(output_dir, paste0("Earth_dicho_PC_std_threshold_", PC_std_threshold, "_", ids, "_", flag, flag_2, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")

  return(list(earth_cont = Earth_cont, earth_dicho = Earth_dicho))
}
