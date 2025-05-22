library("data.table")
library(earth)
library(caret)

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

run_earth <- function(ids,
                      PC_std_threshold,
                      size = 5000,
                      degree = 3,
                      blocks,
                      traits_dir,
                      LDpred2_model = c("auto", "grid"),
                      beta_dir,
                      outcome_db,
                      mask = c("mask", "no_mask"),
                      mask_dir = NULL) {
  if (!file.exists(file.path(blocks))) {
    stop("Blocks file does not exist.")
  }

  if (!file.exists(file.path(cross_val_PC_SD))) {
    stop("Cross validation PC SD file does not exist.")
  }

  # Mar to Earth predict function, no sign
  earth_predict <- function(beta_disc, beta_val, test_data, ind) {
    x <- as.data.frame(beta_disc)

    # 	nk = 5 * dim(x)[2]

    cv_mars_beta <- earth(x = x, y = beta_val)
    # 	cv_mars_beta <- earth(x = x, y = beta_val, thresh = 0.001)
    # 	cv_mars_beta <- earth(x = x, y = beta_val, degree = 3, thresh = 0.001)
    # 	cv_mars_beta <- earth(x = x, y = beta_val, degree = degree, nk = nk, trace = FALSE, thresh = 0.0001)

    # 	save(cv_mars_beta, file = paste0("Model/", colnames(beta_validation)[ind], "_earth_PC_std_threshold_", PC_std_threshold, "_5fold_CV_id_", ids,".RData" ))

    pred <- predict(cv_mars_beta, newdata = test_data)
    return(pred)
  }


  # Main function
  setwd(paste0(trait_dir, "/Traits_", LDpred2_model, "/"))

  block <- as.matrix(fread(blocks))
  # mask_trait_id = as.matrix(fread("Data/masked_traits_20240409.txt", header = F)) #include the GWAS info which will be masked

  beta_validation <- as.numeric()
  beta_discovery <- as.numeric()
  PC_SD <- as.numeric()

  for (chr in 1:22)
  {
    for (set in 1:block[chr])
    {
      # 		ukb_dis_beta = as.matrix(fread(paste0("Betas/20240207/UKB_LDPRED2_PC_Betas_chr_", chr, "_", set, "_disc.txt"), header = T))
      # 		ukb_dis_beta = as.matrix(fread(paste0("Betas/20240207/UKB_Univariate_PC_Betas_chr_", chr, "_", set, "_disc.txt"), header = T))
      ukb_dis_beta <- as.matrix(fread(paste0("Betas/", outcome_db, "_LDPRED2_PC_Betas_chr_", chr, "_", set, "_disc.txt"), header = T))
      PC_SD <- c(PC_SD, ukb_dis_beta[, 1])

      beta_validation <- rbind(beta_validation, ukb_dis_beta[, -1])

      ldpred2_beta <- as.matrix(fread(paste0("Betas/GWAS_LDPRED2_PC_Betas_chr_", chr, "_", set, "_disc.txt"), header = T))
      beta_discovery <- rbind(beta_discovery, ldpred2_beta[, -1])
    }
  }
  dim(beta_discovery)

  PC_SD_Data <- as.matrix(fread(cross_val_PC_SD))
  index_test <- which(PC_SD_Data[, 4] == ids & PC_SD_Data[, 5] >= PC_std_threshold)
  index_train <- which(PC_SD_Data[, 4] != ids & PC_SD_Data[, 5] >= PC_std_threshold)

  beta_validation_train <- beta_validation[index_train, ]

  predict_earth <- as.numeric()
  for (i in 1:dim(beta_validation)[2])
  {
    beta_discovery_train <- cbind(beta_discovery[index_train, ], PC_SD[index_train])
    test_data <- cbind(beta_discovery[index_test, ], PC_SD[index_test])

    if (file.exists(paste0(mask_dir, colnames(beta_validation)[i], "_masked.txt"))) {
      mask_gwas <- as.matrix(fread(paste0(mask_dir, colnames(beta_validation)[i], "_masked.txt"), header = F))
      mask_ids <- match(mask_gwas, colnames(beta_discovery_train))

      beta_discovery_train <- beta_discovery_train[, -mask_ids]
      test_data <- test_data[, -mask_ids]
      # 		print(i)
    }
    pred_earth <- earth_predict(beta_discovery_train, beta_validation_train[, i], test_data, i)

    predict_earth <- cbind(predict_earth, pred_earth)
    print(paste(i, "is done"))
  }
  colnames(predict_earth) <- colnames(beta_validation)

  # Get opt validation PRS
  if (file.exists(paste0("Geno_disc_PCA/Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_val.txt"))) {
    geno_PCA <- as.matrix(fread(paste0("Geno_disc_PCA/Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_val.txt"), header = F))
  } else {
    geno_PCA <- as.matrix(fread(paste0("Geno_disc_PCA/Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", ids, "_val.txt"), header = F))
  }
  earth_PRS <- geno_PCA %*% predict_earth

  rm(geno_PCA)
  gc()

  write.table(earth_PRS, paste0(output_dir, "/Earth_PRS_PC_std_threshold_", PC_std_threshold, "_index_", ids, "_val.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  val_PRs <- earth_PRS

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

  return(list(val_earth_PRS = val_PRS, disc_earth_PRS = disc_PRS))
}
