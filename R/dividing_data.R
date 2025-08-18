#' @param blocks_info File containing three columns of info pertaining to <chr> <set> <id> of the rotation files
#' @param rotations_dir Path to the directory containing rotation files. Each file must be of the form G_PC_SD_chr_<chr>_set_<set>_id_<id>.RData and contain a single object named G_PC_SD.
#' @param Output directory where "Cross_validations_groups.txt" and "Cross_validation_PC_SD.txt" 
#' @return Returns grousp for CV with matching SD
#' @export

library("data.table")
library("stringr")

cross_val_PC_SD <- function(blocks_info,
			    rotation_dir,
			    output_dir){
blocks_orig = as.matrix(fread(blocks_info, header = F))

#blocks1 = matrix(unlist(strsplit(blocks_orig[, 1], "_")), nrow=length(blocks), byrow = TRUE)
#blocks1[,8] = str_replace(blocks1[,8], ".RData", "")
#blocks = blocks1[, c(4, 6, 8)]
blocks <- blocks_orig
colnames(blocks) = c("chr", "set", "ids")
class(blocks) = "numeric"

set.seed(100)

random_vec <- sample(1:5, dim(blocks)[1], prob = c(0.2, 0.2, 0.2, 0.2, 0.2), replace=TRUE)
blocks = cbind(blocks, random_vec)
colnames(blocks)[4] = "groups"

blocks_sort = blocks[order(blocks[,1], blocks[,2]), ]

write.table(blocks_sort, paste0(output_dir,"/Cross_validation_groups.txt"), col.names=T, row.names=F, quote=F, sep="\t")		

PC_SD_Data = as.numeric()
for(i in 1:dim(blocks_sort)[1])
{
	item = blocks_sort[i,]
	load(paste0(rotation_dir,"/G_PC_SD_chr_", item[1], "_set_", item[2], "_id_", item[3], ".RData"))
	temp = cbind(item[1], item[2], item[3], item[4], G_PC_SD)
	PC_SD_Data = rbind(PC_SD_Data, temp)
}

colnames(PC_SD_Data)[1:4] = colnames(blocks_sort)
write.table(PC_SD_Data, paste0(output_dir,"/Cross_validation_PC_SD.txt"), col.names=T, row.names=F, quote=F, sep="\t")		

}
