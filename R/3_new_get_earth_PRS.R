#' @param ids Fold id number (1 to 5 for 5-fold cross validation.)
#' @param PC_std_threshold PC SD threshold for filtering of rotated genotype data
#' @param size Number of continguous variants per block. Default set to 5000.
#' @param degree Degree parameter set for earth to determine allowed interactions. Default set to 3.
#' @param blocks Path to the blocks file.
#' @param traits_dir Path to the directory where the LDpred2 output was saved.
#' @param LDpred2_model Use "auto" (no phenotype data available) or "grid" (phenotype data available) for LDpred2.
#' @param outcome_db Name of outcome database (e.g. UKB).
#' @param mask Specifies whether or not to applying masking condition (leaving out most highly associated trait)
#' @param mask_dir Path to directory containing data for masking. Not required if mask set to "no_mask". Files must be in the format *_masked.txt
#' @param output_dir Path to output directory for earth PRS.
#' @return Returns a list of two elements: disc_earth_PRS and val_earth_PRS
#' @export
run_earth <- function(ids,
                      PC_std_threshold,
                      size = 5000,
                      degree = 3,
                      blocks,
                      cross_val_PC_SD,
                      traits_dir,
                      LDpred2_model,
                      beta_dir,
                      outcome_db,
                      mask = "no_mask",
                      mask_dir = NULL,
                      output_dir) {
  suppressMessages(library("data.table"))
  suppressMessages(library(earth))
  suppressMessages(library(caret))

  blocks <- file.path(blocks)
  cross_val_PC_SD <- file.path(cross_val_PC_SD)
  if (!file.exists(blocks)) {
    stop("Blocks file does not exist.")
  }
  if (!file.exists(cross_val_PC_SD)) {
    stop("Cross validation PC SD file does not exist.")
  }
  if (LDpred2_model != "auto" && LDpred2_model != "grid") {
    stop("LDpred2_model must be either 'auto' or 'grid'.")
  }
  if (mask != "mask" && mask != "no_mask") {
    stop("Mask must be either 'mask' or 'no_mask'.")
  }

  # Mar to Earth predict function, no sign
  earth_predict <- function(beta_disc, beta_val, test_data, ind) {
    x <- as.data.frame(beta_disc)

    cv_mars_beta <- earth(x = x, y = beta_val)

    pred <- predict(cv_mars_beta, newdata = test_data)
    return(pred)
  }

  # Main function
  dir <- file.path(traits_dir, paste0("Traits_", LDpred2_model))

  block <- as.matrix(fread(blocks))

  beta_validation <- as.numeric()
  beta_discovery <- as.numeric()
  PC_SD <- as.numeric()

  for (chr in 1:22)
  {
    for (set in 1:block[chr])
    {
      ukb_dis_beta <- as.matrix(fread(file.path(dir, "Betas", paste0(outcome_db, "_LDPRED2_PC_Betas_chr_", chr, "_", set, "_disc.txt")), header = T))
      PC_SD <- c(PC_SD, ukb_dis_beta[, 1])

      beta_validation <- rbind(beta_validation, ukb_dis_beta[, -1])

      ldpred2_beta <- as.matrix(fread(file.path(dir, "Betas", paste0("GWAS_LDPRED2_PC_Betas_chr_", chr, "_", set, "_disc.txt")), header = T))
      beta_discovery <- rbind(beta_discovery, ldpred2_beta[, -1])
    }
  }
  dim(beta_discovery)

  PC_SD_Data <- as.matrix(fread(cross_val_PC_SD))
  index_test <- which(PC_SD_Data[, 4] == ids & PC_SD_Data[, 5] >= PC_std_threshold)
  index_train <- which(PC_SD_Data[, 4] != ids & PC_SD_Data[, 5] >= PC_std_threshold)

  beta_validation_train <- beta_validation[index_train, ]

  predict_earth <- as.numeric()
  for (i in seq_len(dim(beta_validation)[2]))
  {
    beta_discovery_train <- cbind(beta_discovery[index_train, ], PC_SD[index_train])
    test_data <- cbind(beta_discovery[index_test, ], PC_SD[index_test])

    if (mask == "mask" && file.exists(file.path(mask_dir, paste0(colnames(beta_validation)[i], "_masked.txt")))) {
      mask_gwas <- as.matrix(fread(file.path(mask_dir, paste0(colnames(beta_validation)[i], "_masked.txt")), header = F))
      mask_ids <- match(mask_gwas, colnames(beta_discovery_train))

      beta_discovery_train <- beta_discovery_train[, -mask_ids]
      test_data <- test_data[, -mask_ids]
    }
    pred_earth <- earth_predict(beta_discovery_train, beta_validation_train[, i], test_data, i)

    predict_earth <- cbind(predict_earth, pred_earth)
    print(paste(i, "is done"))
  }
  colnames(predict_earth) <- colnames(beta_validation)

  # Get opt validation PRS
  if (file.exists(file.path(dir, "Geno_disc_PCA", paste0("Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_val.txt")))) {
    geno_PCA <- as.matrix(fread(file.path(dir, "Geno_disc_PCA", paste0("Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_val.txt")), header = F))
  } else {
    geno_PCA <- as.matrix(fread(file.path(dir, "Geno_disc_PCA", paste0("Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_val.txt")), header = F))
  }
  earth_PRS <- geno_PCA %*% predict_earth

  rm(geno_PCA)
  gc()

  write.table(earth_PRS, paste0(output_dir, "/Earth_PRS_PC_std_threshold_", PC_std_threshold, "_index_", ids, "_val.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  val_PRS <- earth_PRS

  # Get true validation PRS
  if (file.exists(paste0("Geno_disc_PCA/Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_disc.txt"))) {
    geno_PCA <- as.matrix(fread(paste0("Geno_disc_PCA/Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_disc.txt"), header = F))
  } else {
    geno_PCA <- as.matrix(fread(paste0("Geno_disc_PCA/Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_disc.txt"), header = F))
  }


  earth_PRS <- geno_PCA %*% predict_earth

  rm(geno_PCA)
  gc()

  write.table(earth_PRS, paste0(output_dir, "/Earth_PRS_PC_std_threshold_", PC_std_threshold, "_index_", ids, "_disc.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  return(list(val_earth_PRS = val_PRS, disc_earth_PRS = earth_PRS))
}
