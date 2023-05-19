## Assessing our SNPs for ones that are incorrect and should be thrown out
## Human genome
human_snps = read.table("/Users/rachelagoglia/Desktop/GRCH38_SNPreview.ALL.SPLIT_SPECIES.txt", header=FALSE)

expLvl_human = human_snps$V3+human_snps$V4+human_snps$V5
expLvl_chimp = human_snps$V6+human_snps$V7+human_snps$V8

percent_human = (human_snps$V3 / expLvl_human) * 100
percent_chimp = (human_snps$V6 / expLvl_chimp) * 100

final = human_snps[which(expLvl_chimp>2 & expLvl_human >2 & percent_human>90 & percent_chimp>90),]

write.table(file="/Users/rachelagoglia/Desktop/GRCh38_final_SNPs.SPLIT_SPECIES.txt", x=final[,1:2], quote=FALSE, col.names = FALSE, row.names = FALSE)

## Chimp genome
chimp_snps = read.table("/Users/rachelagoglia/Desktop/PanTro5_SNPreview.ALL.SPLIT_SPECIES.ALLCHR.txt", header=FALSE)

expLvl_human = chimp_snps$V3+chimp_snps$V4+chimp_snps$V5
expLvl_chimp = chimp_snps$V6+chimp_snps$V7+chimp_snps$V8

percent_human = (chimp_snps$V3 / expLvl_human) * 100
percent_chimp = (chimp_snps$V6 / expLvl_chimp) * 100

final = chimp_snps[which(expLvl_chimp>2 & expLvl_human >2 & percent_human>90 & percent_chimp>90),]

write.table(file="/Users/rachelagoglia/Desktop/PanTro5_final_SNPs.SPLIT_SPECIES.txt", x=final[,1:2], quote=FALSE, col.names = FALSE, row.names = FALSE)

