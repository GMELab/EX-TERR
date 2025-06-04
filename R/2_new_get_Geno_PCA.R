library("data.table")

#' @param id Fold id number, from 1 to 5
#' @param flag Type of set: disc (discovery) or val (validation)
#' @param PC_std_threshold PC SD threshold for filtering of rotated genotype data
#' @param genotype_dir Directory containing genotype directions (either Geno_disc or Geno_val)
#' @param outcome_db Name of outcome database (e.g. UKB)
#' @param rotations_dir Path to directory containing rotations G_PC_SD_chr_<chr>_set_<set>_id_<group_id>.RData.
#' @param cv_groups Path to cross validation groups file (Cross_validation_groups.txt)
#' @param output_dir Name of target output dir
#' @param return Element geno_pc.


geno_PCA <- function(id,
                     flag = c("disc", "val"),
                     PC_std_threshold,
                     genotype_dir,
                     outcome_db,
                     rotations_dir,
                     cv_groups,
                     output_dir = NULL) {
  if (!file.exists(file.path(cv_groups))) {
    stop("Cross validation groups file does not exist.")
  }

  get_PCA_geno <- function(item) {
    chr <- item[1]
    set <- item[2]
    group_id <- item[3]

    load(file.path(genotype_dir, paste0("Geno_", flag, "/UKB_09_", chr, "_", set, ".RData"))) # load genotypes


    if (grepl(flag, rotations_dir)) { # Ensure correct set: disc or val
      load(file.path(rotations_dir, paste0(outcome_db, "_09_", chr, "_", set, ".RData"))) # load genotypes
    }

    starting <- (group_id - 1) * 5000 + 1
    ending <- group_id * 5000
    if (ending > dim(geno_data)[2]) {
      ending <- dim(geno_data)[2]
    }

    G <- geno_data[, starting:ending]

    # load G_PC_rotation. The number of PCs is 5000 in one block
    load(file.path(rotations_dir, paste0("G_PC_chr_", chr, "_set_", set, "_id_", group_id, ".RData")))
    # load G_PC_SD
    load(file.path(rotations_dir, paste0("G_PC_SD_chr_", chr, "_set_", set, "_id_", group_id, ".RData")))

    select_PC_1 <- which(G_PC_SD >= PC_std_threshold)
    G_PC_rotation <- G_PC_rotation[, select_PC_1]

    Final_G_PCs <- G %*% G_PC_rotation

    geno_data_PCA <- scale(Final_G_PCs, center = T, scale = G_PC_SD[select_PC_1])

    rm(G)
    rm(G_PC_rotation)
    rm(Final_G_PCs)
    rm(geno_data)
    gc()

    return(geno_data_PCA)
  }


  # Main function

  blocks_orig <- as.matrix(fread(cv_groups))
  blocks <- blocks_orig[which(blocks_orig[, 4] == id), ]

  geno_PCA <- as.numeric()
  for (i in seq_len(dim(blocks)[1]))
  {
    item <- blocks[i, ]

    geno_data_PCA <- get_PCA_geno(item)
    geno_PCA <- cbind(geno_PCA, geno_data_PCA)

    rm(geno_data_PCA)
    gc()

    print(paste("ids", id, "number of ", i, "is done"))
  }

  if (!is.null(output_dir)) {
    write.table(geno_PCA, file.path(output_dir, paste0("Geno_PCA_PC_std_threshold_", PC_std_threshold, "_", id, "_", flag, ".txt")), col.names = F, row.names = F, quote = F, sep = "\t")
  }

  return(geno_pc = geno_PCA)
}
