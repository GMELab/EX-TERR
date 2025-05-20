library("data.table")

#' @param chr Chromsome number
#' @param size Number of continguous variants per block. Default set to 5000.
#' @param rotations_dir Path to the directory containing rotation files. Each file must be of the form G_PC_SD_chr_<chr>_set_<set>_id_<id>.RData and contain a single object named G_PC_SD.
#' @param blocks Path to the blocks file.
#' @param outcome_db Name of outcome database (e.g. UKB).
#' @param genotype_dir Path to directory containing genotype .Rdata files. Each file must be of the form <outcome_db>_09_<chr>_<set>.RData and contain a single object named Data.
#' @return A list containing two elements: g_pc_rotation  and g_pc_sd.


get_rotations < function(chr,
                         size = 5000,
                         rotations_dir,
                         blocks,
                         genotype_dir){

  if (!file.exists(file.path(blocks))) {
    stop("Blocks file does not exist.")
  }

get_PCA_geno = function(Data)
{
  cycles = ceiling(dim(Data)[2] / size)
  ending = 0
  block_end = 0

  for(i in 1:cycles)
  {
    starting = ending + 1
    ending = ending + size
    if(ending > dim(Data)[2])
    {
      ending = dim(Data)[2]
    }

    if(!file.exists(paste0(rotations_dir,"G_PC_chr_", chr, "_set_", set, "_id_", i, ".RData")))
    {
      G = Data[, starting:ending]
      G_PC <- prcomp( G , scale = T )

      G_PC_rotation <- G_PC$rotation
      G_PC_SD <- G_PC$sdev
      save(G_PC_rotation, file = paste0(rotations_dir,"G_PC_chr_", chr, "_set_", set, "_id_", i, ".RData"))
      save(G_PC_SD, file = paste0(rotations_dir,"G_PC_SD_chr_", chr, "_set_", set, "_id_", i, ".RData"))

      rm(G)
      rm(G_PC)
      gc()
    }
  }
  return(1)
}


# Main function
#setwd("/genetics3/maos/Geno_PC_external_GWAS")
block = as.matrix(fread(blocks))

for (set in 1:block[chr])
 {
  load(paste0(genotype_dir, "/", outcome_db, "_09_", chr, "_", set, ".RData")) #load genotypes

  geno_data_PCA = get_PCA_geno(geno_data)
  rm(geno_data)
  gc()
 }
 return(g_pc_rotation = G_PC_rotation, g_pc_sd = G_PC_SD)
}

