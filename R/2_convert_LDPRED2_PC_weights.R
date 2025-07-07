#' @param trait_type Specifies the type of trait: auto, grid or outcome
#' @param chr chromosome number
#' @param flag Type of set: disc (discovery) or val (validation)
#' @param size Size/number of rotated components per block. Default is 5000.
#' @param outcome_db Name of outcome database (e.g. UKB).
#' @param rotations_dir Path to the directory containing rotation files. Each file must be of the form G_PC_SD_chr_<chr>_set_<set>_id_<id>.RData and contain a single object named G_PC_SD.
#' @param LDpred2_model Use "auto" (no phenotype data available) or "grid" (phenotype data available) for LDpred2.
#' @param traits_dir Path to the directory containing LDpred2 output sorted by traits.
#' @param trait_list_dir Path to directory that contains trait list files. Must be in the format "GWas_list_<LDpred2_model>"contain file specifying gwas_list_auto.txt and gwas_list_grid.txt.
#' @param blocks Path to the blocks file.
#' @param bim_dir Path to directory containing outcome genotype .bim files divided into chromsome and set. Must be in the form <outcome_db>_09_<chr>_<set>.bim
#' @return Returns element beta
#' @export
convert_LDpred2 <- function(trait_type,
                            chr,
                            flag,
                            size = 5000,
                            outcome_db,
                            rotations_dir,
                            LDpred2_model,
                            traits_dir,
                            trait_list_dir,
                            blocks,
                            bim_dir,
                            bim_file) {
  suppressMessages(library("data.table"))

  if (!file.exists(file.path(blocks))) {
    stop("Blocks file does not exist.")
  }
  if (!file.exists(file.path(bim_dir))) {
    stop("bim_dir does not exist.")
  }

  get_PC_weights <- function(adj_betas) {
    PC_beta_weight <- as.numeric()

    cycles <- ceiling(dim(adj_betas)[1] / size)
    ending <- 0

    for (i in 1:cycles)
    {
      starting <- ending + 1
      ending <- ending + size
      if (ending > dim(adj_betas)[1]) {
        ending <- dim(adj_betas)[1]
      }

      # load G_PC_rotation. The number of PCs is 5000 in one block
      load(file.path(rotations_dir, paste0("G_PC_chr_", chr, "_set_", set, "_id_", i, ".RData")))
      # load G_PC_SD
      load(file.path(rotations_dir, paste0("G_PC_SD_chr_", chr, "_set_", set, "_id_", i, ".RData")))
      weight_block <- t(t(adj_betas[starting:ending, ]) %*% G_PC_rotation)
      weight_block_adj <- weight_block * G_PC_SD

      PC_beta_weight <- rbind(PC_beta_weight, cbind(G_PC_SD, weight_block_adj))
    }
    return(PC_beta_weight)
  }

  collect_LDPRED2_betas <- function(traits, tag) {
    GWAS_LDPRED2_betas <- as.numeric()
    exists <- c()
    for (i in seq_along(traits))
    {
      trait <- traits[i]
      path <- file.path(traits_dir, paste0("Traits_", LDpred2_model), trait, "Betas", "LDpred2_betas_signed_adj.txt")
      if (!file.exists(path)) {
        next
      }
      exists <- c(exists, trait)
      gwas <- as.matrix(fread(path))
      GWAS_LDPRED2_betas <- cbind(GWAS_LDPRED2_betas, gwas[, 3])
    }
    GWAS_LDPRED2_betas <- matrix(GWAS_LDPRED2_betas, ncol = length(exists), byrow = FALSE)
    colnames(GWAS_LDPRED2_betas) <- exists
    class(GWAS_LDPRED2_betas) <- "numeric"
    # save(GWAS_LDPRED2_betas, file = paste0("/genetics3/maos/Geno_PC_external_GWAS/Traits_", tag, "/GWAS_LDPRED2_betas.RData"))
    save(GWAS_LDPRED2_betas, file = file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))

    return(GWAS_LDPRED2_betas)
  }


  # Main function
  # setwd("/genetics_work2/lea/10-12-EX-TERR")
  block <- as.matrix(fread(blocks, header = F))

  if (LDpred2_model == "auto") {
    traits_auto <- as.matrix(fread(file.path(trait_list_dir, paste0("Gwas_list_", LDpred2_model, ".txt")), header = F))

    if (!file.exists(file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))) {
      ldpred2_beta <- collect_LDPRED2_betas(traits_auto, "auto")
    } else {
      # load(paste0("/genetics3/maos/Geno_PC_external_GWAS/Traits_auto/GWAS_LDPRED2_betas.RData"))
      load(file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))

      ldpred2_beta <- GWAS_LDPRED2_betas
    }
  } else if (LDpred2_model == "grid") {
    traits_grid <- as.matrix(fread(file.path(trait_list_dir, paste0("Gwas_list_", LDpred2_model, ".txt")), header = F))

    if (!file.exists(file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))) {
      ldpred2_beta <- collect_LDPRED2_betas(traits_grid, "grid")
    } else {
      # load(paste0("/genetics3/maos/Geno_PC_external_GWAS/Traits_GRID/GWAS_LDPRED2_betas.RData"))
      load(file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))
      ldpred2_beta <- GWAS_LDPRED2_betas
    }
  } else if (LDpred2_model == "outcome") {
    traits_outcome <- as.matrix(fread(file.path(trait_list_dir, "outcome_list_final.txt"), header = F))

    if (!file.exists(file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))) {
      ldpred2_beta <- collect_LDPRED2_betas(traits_outcome, outcome_db)
    } else {
      load(file.path(traits_dir, paste0("Traits_", LDpred2_model), "GWAS_LDPRED2_betas.RData"))
      ldpred2_beta <- GWAS_LDPRED2_betas
    }
  }

  bim <- as.matrix(fread(bim_file, header = F))
  print(dim(ldpred2_beta))
  ldpred2_beta_chr <- ldpred2_beta[which(as.numeric(bim[, 1]) == chr), , drop = FALSE]
  print(dim(ldpred2_beta_chr))

  if (!dir.exists(file.path(traits_dir, paste0("Traits_", LDpred2_model), "Betas"))) {
    dir.create(file.path(traits_dir, paste0("Traits_", LDpred2_model), "Betas"), recursive = TRUE)
  }

  ending <- 0
  for (set in 1:block[chr])
  {
    bim <- as.matrix(fread(file.path(bim_dir, paste0("Geno_", flag), paste0(outcome_db, "_09_", chr, "_", set, ".bim")), header = F))
    starting <- ending + 1
    ending <- ending + dim(bim)[1]

    betas <- get_PC_weights(ldpred2_beta_chr[starting:ending, , drop = FALSE])

    if (trait_type == "outcome") {
      write.table(betas, file = file.path(traits_dir, paste0("Traits_", LDpred2_model), "Betas", paste0(outcome_db, "_LDPRED2_PC_Betas_chr_", chr, "_", set, "_", flag, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")
    } else {
      write.table(betas, file = file.path(traits_dir, paste0("Traits_", LDpred2_model), "Betas", paste0("GWAS_LDPRED2_PC_Betas_chr_", chr, "_", set, "_", flag, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")
    }
  }
  return(list(beta = betas))
}
