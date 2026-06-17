#' @param ids Fold id number (1 to 5 for 5-fold cross validation.)
#' @param cv_group Path to file containing cross_validation_groups.txt
#' @param outcome_list Path to file containing outcome list.
#' @param mask_dir Path to directory containing data for masking. Not required if mask set to "no_mask". Files must be in the format *_masked.txt
#' @param pheno_dir Path to directory containing phenotype information separated by trait. Files are named <trait>_pheno.txt
#' @return Returns a list containing four elements: PRS_cont_train, PRS_cont_test, PRS_dicho_train, PRS_dicho_test
#' @export
LASSO_cv <- function(ids,
                     cv_group,
                     outcome_list,
                     mask,
                     mask_dir = NULL,
                     pheno_dir,
                     prs_dir) {
  suppressMessages(library("data.table"))
  suppressMessages(library("glmnet"))

  # Main function

  blocks_orig <- as.matrix(fread(cv_group)) # Randomly divided into 5 groups
  testing_traits <- as.matrix(fread(outcome_list, header = F)) # 71 test traits
  # mask_trait_id = as.matrix(fread("Data/Mask_traits.txt", header = F)) #include the GWAS info which will be masked

  train_blocks <- blocks_orig[which(blocks_orig[, 4] != ids), ]
  test_blocks <- blocks_orig[which(blocks_orig[, 4] == ids), ]

  # Get train_PRS and test_PRS
  train_prs <- as.matrix(fread(file.path(prs_dir, "PRS_chr_1_1_1_disc.txt")))
  train_PRS <- matrix(0, ncol = dim(train_prs)[2], nrow = dim(train_prs)[1])
  for (i in seq_len(dim(train_blocks)[1]))
  {
    item <- train_blocks[i, ]
    train_prs <- as.matrix(fread(file.path(prs_dir, paste0("PRS_chr_", item[1], "_", item[2], "_", item[3], "_disc.txt"))))
    train_PRS <- train_PRS + train_prs
  }

  one_out_prs <- as.matrix(fread(file.path(prs_dir, "PRS_chr_1_1_1_disc.txt"))) # This is one out group in discovery set
  one_out_PRS <- matrix(0, ncol = dim(one_out_prs)[2], nrow = dim(one_out_prs)[1])
  for (i in seq_len(dim(test_blocks)[1]))
  {
    item <- test_blocks[i, ]
    one_out_prs <- as.matrix(fread(file.path(prs_dir, paste0("PRS_chr_", item[1], "_", item[2], "_", item[3], "_disc.txt"))))
    one_out_PRS <- one_out_PRS + one_out_prs
  }

  test_prs <- as.matrix(fread(file.path(prs_dir, "PRS_chr_1_1_1_val.txt"))) # This is one out group in validation set
  test_PRS <- matrix(0, ncol = dim(test_prs)[2], nrow = dim(test_prs)[1])
  for (i in seq_len(dim(test_blocks)[1]))
  {
    item <- test_blocks[i, ]
    test_prs <- as.matrix(fread(file.path(prs_dir, paste0("PRS_chr_", item[1], "_", item[2], "_", item[3], "_val.txt"))))
    test_PRS <- test_PRS + test_prs
  }

  # LASSO training
  continuous_trait <- as.character()
  one_out_PRS_cont_train <- as.numeric()
  one_out_PRS_cont_test <- as.numeric()

  dicho_trait <- as.character()
  one_out_PRS_dicho_train <- as.numeric()
  one_out_PRS_dicho_test <- as.numeric()

  for (i in seq_len(length(testing_traits)))
  {
    trait <- testing_traits[i]

    x_train <- train_PRS # in disc set
    one_train <- one_out_PRS # one out in disc set
    one_test <- test_PRS # one out in val set

    if (mask == "mask" && file.exists(file.path(mask_dir, paste0(trait, "_masked.txt")))) {
      mask_gwas <- as.matrix(fread(file.path(mask_dir, paste0(trait, "_masked.txt")), header = F))
      mask_ids <- match(mask_gwas, colnames(train_PRS))

      x_train <- x_train[, -mask_ids, drop = FALSE]
      one_train <- one_train[, -mask_ids, drop = FALSE]
      one_test <- one_test[, -mask_ids, drop = FALSE]
    }

    # phenos = as.matrix(fread(paste0("/genetics3/maos/Geno_PC_external_GWAS/Traits_UKB/", trait, "/Pheno/", trait, "_disc_updated.txt"), header = T))
    phenos <- as.matrix(fread(file.path(pheno_dir, paste0(trait, "_disc.txt")), header = T))
    # define response variable
    y_train <- phenos[, 3, drop = F]
    if (max(y_train) != 1 || min(y_train) != 0) {
      if (ncol(x_train) > 1) {
        # perform k-fold cross-validation to find optimal lambda value
        cv_model <- cv.glmnet(x_train, y_train, alpha = 1)
        # find optimal lambda value that minimizes test MSE
        best_lambda <- cv_model$lambda.min
        best_model <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)

        lasso_predict_train <- predict(best_model, s = 0.2, newx = one_train)
        lasso_predict_test <- predict(best_model, s = 0.2, newx = one_test)
      } else {
        best_model <- glm(y_train ~ x_train, family = "gaussian")
        lasso_predict_train <- predict(best_model, newdata = data.frame(x_train = one_train))
        lasso_predict_test <- predict(best_model, newdata = data.frame(x_train = one_test))
      }

      one_out_PRS_cont_train <- cbind(one_out_PRS_cont_train, lasso_predict_train)
      one_out_PRS_cont_test <- cbind(one_out_PRS_cont_test, lasso_predict_test)

      continuous_trait <- c(continuous_trait, trait)
    } else {
      if (ncol(x_train) > 1) {
        # perform k-fold cross-validation to find optimal lambda value
        cv_model <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial")

        # find optimal lambda value that minimizes test MSE
        best_lambda <- cv_model$lambda.min
        best_model <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda, family = "binomial")

        lasso_predict_train <- predict(best_model, s = 0.2, newx = one_train)
        lasso_predict_test <- predict(best_model, s = 0.2, newx = one_test)
      } else {
        best_model <- glm(y_train ~ x_train, family = "binomial")
        lasso_predict_train <- predict(best_model, newdata = data.frame(x_train = one_train))
        lasso_predict_test <- predict(best_model, newdata = data.frame(x_train = one_test))
      }

      one_out_PRS_dicho_train <- cbind(one_out_PRS_dicho_train, lasso_predict_train)
      one_out_PRS_dicho_test <- cbind(one_out_PRS_dicho_test, lasso_predict_test)

      dicho_trait <- c(dicho_trait, trait)
    }
    print(paste("IDs", ids, i, trait, "is done"))
  }

  one_out_PRS_cont_train <- as.matrix(one_out_PRS_cont_train)
  one_out_PRS_cont_test <- as.matrix(one_out_PRS_cont_test)
  one_out_PRS_dicho_train <- as.matrix(one_out_PRS_dicho_train)
  one_out_PRS_dicho_test <- as.matrix(one_out_PRS_dicho_test)
  print(continuous_trait)

  colnames(one_out_PRS_cont_train) <- continuous_trait
  colnames(one_out_PRS_cont_test) <- continuous_trait

  colnames(one_out_PRS_dicho_train) <- dicho_trait
  colnames(one_out_PRS_dicho_test) <- dicho_trait


  write.table(one_out_PRS_cont_train, file.path(prs_dir, paste0("R6_LASSO_continuous_disc_", ids, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(one_out_PRS_cont_test, file.path(prs_dir, paste0("R6_LASSO_continuous_val_", ids, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")

  write.table(one_out_PRS_dicho_train, file.path(prs_dir, paste0("R6_LASSO_dicho_disc_", ids, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(one_out_PRS_dicho_test, file.path(prs_dir, paste0("R6_LASSO_dicho_val_", ids, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")

  return(list(PRS_cont_train = one_out_PRS_cont_train, PRS_cont_test = one_out_PRS_cont_test, PRS_dicho_train = one_out_PRS_dicho_train, PRS_dicho_test = one_out_PRS_dicho_test))
}
