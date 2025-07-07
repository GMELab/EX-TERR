#' @param chr chromosome number
#' @param flag Type of set: disc (discovery) or val (validation)
#' @param size Size/number of rotated components per block. Default is 5000.
#' @param traits_dir Path to the directory containing LDpred2 output sorted by traits.
#' @param genotype_dir Path to directory containing genotype files. Directory contains directories Geno_<flag> and files <outcome_db>_09_<chr>_<set>.Rdata
#' @param blocks Path to the blocks file.
#' @param trait_list_dir Path to directory containing traits_list sorted by ldpred_model: gwas_list_auto.txt and gwas_list_grid.txt
#' @export
get_block_PRS <- function(chr,
                          flag,
                          size = 5000,
                          traits_dir,
                          genotype_dir,
                          blocks,
                          trait_list_dir,
                          outcome_db,
                          output_dir,
                          LDpred2_model) {
  suppressMessages(library("data.table"))

  if (flag != "disc" && flag != "val") {
    stop("Flag must be either 'disc' or 'val'.")
  }
  if (!file.exists(file.path(blocks))) {
    stop("Blocks file does not exist.")
  }


  collect_LDPRED2_betas <- function(traits, tag) {
    GWAS_LDPRED2_betas <- as.numeric()
    for (i in seq_along(traits))
    {
      trait <- traits[i]
      gwas <- as.matrix(fread(file.path(traits_dir, paste0("Traits_", tag), trait, "Betas", "LDpred2_betas_signed_adj.txt")))
      GWAS_LDPRED2_betas <- cbind(GWAS_LDPRED2_betas, gwas[, 3])
    }
    colnames(GWAS_LDPRED2_betas) <- traits
    class(GWAS_LDPRED2_betas) <- "numeric"

    return(GWAS_LDPRED2_betas)
  }


  get_PRS <- function(adj_betas) {
    cycles <- ceiling(dim(adj_betas)[1] / size)
    ending <- 0
    load(file.path(genotype_dir, paste0("Geno_", flag), paste0(outcome_db, "_09_", chr, "_", set, ".RData")))

    for (i in 1:cycles)
    {
      starting <- ending + 1
      ending <- ending + size
      if (ending > dim(adj_betas)[1]) {
        ending <- dim(adj_betas)[1]
      }
      PRS <- geno_data[, starting:ending] %*% adj_betas[starting:ending, ]

      write.table(PRS, file = file.path(output_dir, paste0("PRS_chr_", chr, "_", set, "_", i, "_", flag, ".txt")), col.names = T, row.names = F, quote = F, sep = "\t")

      print(paste("chr", chr, i, "out of", cycles, " beta is done"))
    }

    rm(geno_data)
    gc()

    return(0)
  }



  # Main function
  block <- as.matrix(blocks, header = F)

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
  } else {
    stop("LDpred2_model must be either 'auto' or 'grid'.")
  }

  bim <- as.matrix(fread(file.path(genotype_dir, "Geno_disc", paste0(outcome_db, "_final.bim")), header = F))

  ldpred2_beta_chr <- ldpred2_beta[which(as.numeric(bim[, 1]) == chr), ]

  SNP_ending <- 0
  for (set in 1:block[chr])
  {
    bim <- as.matrix(fread(file.path(genotype_dir, paste0("Geno_", flag), paste0(outcome_db, "_09_", chr, "_", set, ".bim")), header = F))
    SNP_starting <- SNP_ending + 1
    SNP_ending <- SNP_ending + dim(bim)[1]

    dummy <- get_PRS(ldpred2_beta_chr[SNP_starting:SNP_ending, ])

    print(paste("chr", chr, "set", set, "is done"))
  }
}
