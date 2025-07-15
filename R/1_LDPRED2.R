#' @param blocks Path to the blocks file.
#' @param trait Name of trait
#' @param outcome_db Name of outcome database (e.g. UKB).
#' @param rds_in .rds file for target outcome genotypes
#' @param bed_in .bed file for target outcome genotypes. Required if .rds file not available.
#' @param sumstats_dir Path to directory containing summary statistics. File should have the name and format of [traits]_sumstats.txt
#' @param freq_dir Path to directory containing MAF files. Each file must be in the form of UKB_09_freq_<chr>_<set>.frq.
#' @param ref_allele_dir Path to the directory containing ref_allele files. Files will be in the form ref_allele_09 and ref_allele_09_from_<outcome>.
#' @param traits_dir Path to the directory where the output files will be saved. If NULL, the files will not be saved.
#' @param trait_type Whether trait is "gwas" or "outcome".
#' @param LDpred2_model Use "auto" (no phenotype data available) or "grid" (phenotype data available) for LDpred2.
#' @param phenotype Path to the phenotype file, required only when using the "grid" model.
#' @param corrections_dir Path to directory containing corrections for Age, Sex and PCs. Files should be in the form <correction>_disc.txt
#' @param ncores Number of cores used.
#' @return A list containing two elements: ldpred2_betadj & ldpred2_beta
#' @export
run_LDPred2 <- function(blocks,
                        trait,
                        outcome_db,
                        rds_in,
                        bed_in = NULL,
                        sumstats_dir,
                        freq_dir,
                        ref_allele_dir,
                        traits_dir = NULL,
                        trait_type,
                        LDpred2_model,
                        phenotype = NULL,
                        corrections_dir = NULL,
                        ncores) {
  suppressMessages(library("data.table"))
  suppressMessages(library(bigsnpr))
  suppressMessages(options(bigstatsr.check.parallel.blas = FALSE))
  suppressMessages(options(default.nproc.blas = NULL))
  suppressMessages(library(magrittr))
  suppressMessages(library(dplyr))
  suppressMessages(library(pscl))
  suppressMessages(library(fmsb))

  phenotype <- file.path(phenotype)
  rds_in <- file.path(rds_in)
  bed_in <- file.path(bed_in)
  blocks <- file.path(blocks)
  if (!file.exists(blocks)) {
    stop("Blocks file does not exist.")
  }
  if (!file.exists(rds_in) && !file.exists(bed_in)) {
    stop("Neither the .rds nor .bed for the outcome genotype data exist. Check path.")
  }
  if (!file.exists(file.path(ref_allele_dir, paste0("ref_allele_09_from_", outcome_db))) || !file.exists(file.path(ref_allele_dir, "ref_allele_09"))) {
    stop("At least one of the ref allele files is missing.")
  }
  if (trait_type != "gwas" && trait_type != "outcome") {
    stop("Trait type must be either 'gwas' or 'outcome'.")
  }
  if (LDpred2_model != "auto" && LDpred2_model != "grid") {
    stop("LDpred2 model must be either 'auto' or 'grid'.")
  }

  get_adj_lassosum_betas <- function(betas) {
    block <- as.matrix(fread(blocks, header = F))

    MAF <- as.numeric()
    for (chr in 1:22)
    {
      for (set in 1:block[chr])
      {
        maf <- as.matrix(fread(file.path(freq_dir, paste0(outcome_db, "_09_freq_", chr, "_", set, ".frq"))))
        MAF <- rbind(MAF, maf)
      }
    }

    length(which(betas[, 2] != MAF[, 2]))
    adj_betas <- as.numeric(betas[, 3]) * sqrt(2 * as.numeric(MAF[, 5]) * (1 - as.numeric(MAF[, 5])))

    betas[, 3] <- adj_betas

    return(betas)
  }


  # Main function
  dir <- file.path(traits_dir, paste0("Traits_", LDpred2_model))
  if (trait_type == "gwas") {
    dir <- file.path(traits_dir, paste0("Traits_", LDpred2_model))
  } else if (trait_type == "outcome") {
    dir <- file.path(traits_dir, paste0("Traits_", outcome_db))
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  if (!dir.exists(file.path(dir, trait, "Betas"))) {
    dir.create(file.path(dir, trait, "Betas"), recursive = TRUE)
  }

  # Grid model
  if (LDpred2_model == "grid") {
    if (!file.exists(file.path(phenotype))) {
      stop("Phenotype information does not exist for Ldpred2 grid model.")
    }

    # Phenotype info
    phenotype <- as.matrix(fread(file.path(phenotype), header = T))
    Age <- as.matrix(fread(file.path(corrections_dir, "Age_disc.txt"), header = T))
    Sex <- as.matrix(fread(file.path(corrections_dir, "Sex_disc.txt"), header = T))
    PCs <- as.matrix(fread(file.path(corrections_dir, "PCs_disc.txt"), header = T))
  }
  sumstats <- bigreadr::fread2(file.path(sumstats_dir, paste0(trait, "_sumstats.txt")))
  # LDpred 2 require the header to follow the exact naming
  names(sumstats) <- c("rsid", "chr", "pos", "a0", "a1", "beta", "beta_se", "N", "p")

  if (!file.exists(rds_in)) {
    # preprocess the bed file (only need to do once for each data set)
    snp_readBed(bed_in)
  }

  # now attach the genotype object
  # 1obj.bigSNP <- snp_attach(paste0("/genetics3/maos/Geno_PC_external_GWAS/Geno_disc/UKB_final.rds"))
  obj.bigSNP <- snp_attach(paste0(rds_in))

  G <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  y <- obj.bigSNP$fam$affection #<- phenotype[, 3]
  class(y) <- "numeric"

  if (LDpred2_model == "auto") {
    sumstats$N[which(is.na(sumstats$N))] <- 30000

    index <- which((sumstats$N) < 0)
    if (length(index) > 0) {
      sumstats$N[index] <- mean(sumstats$N[-index])
    }
  }

  sumstats$n_eff <- sumstats$N
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
  sumstats$chr <- as.numeric(sumstats$chr)
  df_beta <- snp_match(sumstats, map) # use rsid instead of pos

  POS2 <- snp_asGeneticPos(CHR, POS, dir = file.path(dir, trait))

  tmp <- tempfile(tmpdir = file.path(dir, trait, "tmp-data"))

  if (!file.exists(file.path(dir, "corr.RData")) || !file.exists(file.path(dir, "ld.RData"))) {
    for (chr in 1:22) {
      ## indices in 'df_beta'
      ind.chr <- which(df_beta$chr == chr)
      ## indices in 'G'
      ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]

      corr0 <- snp_cor(G,
        ind.col = ind.chr2, size = 3 / 1000,
        infos.pos = POS2[ind.chr2], ncores = ncores
      )

      if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
      print(chr)
    }
    save(corr, file = file.path(dir, "corr.RData"))
    save(ld, file = file.path(dir, "ld.RData"))
  } else {
    load(file.path(dir, "corr.RData"))
    load(file.path(dir, "ld.RData"))
  }

  ldsc <- with(df_beta, snp_ldsc(ld, length(ld),
    chi2 = (beta / beta_se)^2,
    sample_size = n_eff, blocks = NULL
  ))

  # LDpred2

  if (LDpred2_model == "auto") {
    # LDpred2-auto: automatic mode
    ldsc_h2_est <- ldsc[["h2"]]

    if (ldsc_h2_est <= 0) {
      ldsc_h2_est <- 0.00001
    }
    coef_shrink <- 0.95 # reduce this up to 0.4 if you have some (large) mismatch with the LD ref

    set.seed(1) # to get the same result every time

    # takes less than 2 min with 4 cores
    multi_auto <- snp_ldpred2_auto(
      corr, df_beta,
      h2_init = ldsc_h2_est,
      vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = ncores,
      #  use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
      allow_jump_sign = FALSE, shrink_corr = coef_shrink
    )

    # str(multi_auto, max.level = 1)
    # str(multi_auto[[1]], max.level = 1)

    range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))

    keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

    beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    LDpred2_betas <- cbind(df_beta$chr, df_beta$rsid, beta_auto)
  } else if (LDpred2_model == "grid") {
    h2_est <- ldsc[["h2"]]

    h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
    h2_seq[which(h2_seq == 0)] <- 1e-6

    p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)

    params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

    # takes less than 2 min with 4 cores
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = ncores)

    pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])

    params$score <- apply(pred_grid, 2, function(x) {
      if (all(is.na(x))) {
        return(NA)
      }

      if (max(y, na.rm = T) == 1) {
        summary(glm(y ~ x + Age[, 2] + Sex[, 2] + PCs[, 2:11], family = "binomial"))$coef["x", 3]
      } else {
        summary(lm(y ~ x + Age[, 2] + Sex[, 2] + PCs[, 2:11]))$coef["x", 3]
      }
    })

    best_beta_grid <- params %>%
      mutate(id = row_number()) %>%
      # filter(sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      print() %>%
      pull(id) %>%
      beta_grid[, .]

    LDpred2_betas <- cbind(df_beta$chr, df_beta$rsid, best_beta_grid)
  } else {
    stop("LDpred2 model not specified. Please specify auto or grid.")
  }


  colnames(LDpred2_betas) <- c("chr", "rsid", "betas")
  write.table(LDpred2_betas, file.path(dir, trait, "Betas", "LDpred2_betas.txt"), col.names = T, row.names = F, quote = F, sep = "\t")


  # sign the betas whose effect allele is different from that in UKB
  UKB_ref <- as.matrix(fread(file.path(ref_allele_dir, paste0("ref_allele_09_from_", outcome_db)), header = F))
  GWAS_ref <- as.matrix(fread(file.path(ref_allele_dir, "ref_allele_09"), header = F))

  index <- (which(UKB_ref[, 2] != GWAS_ref[, 2]))

  if (length(which(UKB_ref[, 1] != LDpred2_betas[, 2])) > 0) {
    print(paste("Something wrong at ", trait, "LDpred betas"))
  }

  LDpred2_betas[index, 3] <- 0 - as.numeric(LDpred2_betas[index, 3])
  write.table(LDpred2_betas, file.path(dir, trait, "Betas", "LDpred2_betas_signed.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  # Standardize the betas
  ldpred2_beta_adj <- get_adj_lassosum_betas(LDpred2_betas)

  write.table(ldpred2_beta_adj, file.path(dir, trait, "Betas", "LDpred2_betas_signed_adj.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

  return(list(ldpred2_betadj = ldpred2_beta_adj, ldpred2_beta = LDpred2_betas))
}
