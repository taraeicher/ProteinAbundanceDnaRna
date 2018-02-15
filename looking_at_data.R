#Load the Dream Challenge data.
breast_dna = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_breast_CNA_sort_common_gene_16884.txt", row.names = 1)
breast_dna_num = breast_dna
indx <- sapply(breast_dna, is.factor)
breast_dna_num[indx] = lapply(breast_dna[indx], function(x) as.numeric(as.character(x)))

breast_proteome = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_breast_proteome_sort_common_gene_10005.txt", row.names = 1, fill = TRUE)
breast_proteome_num = breast_proteome
indx <- sapply(breast_proteome, is.factor)
breast_proteome_num[indx] = lapply(breast_proteome[indx], function(x) as.numeric(as.character(x)))

breast_rna = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_breast_RNA_sort_common_gene_15107.txt", row.names = 1, fill = TRUE)
breast_rna_num = breast_rna
indx <- sapply(breast_rna, is.factor)
breast_rna_num[indx] = lapply(breast_rna[indx], function(x) as.numeric(as.character(x)))

ova_microarray = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_ova_array_sort_common_gene_15121.txt", row.names = 1)
ova_microarray_num = ova_microarray
indx <- sapply(ova_microarray, is.factor)
ova_microarray_num[indx] = lapply(ova_microarray[indx], function(x) as.numeric(as.character(x)))

ova_dna = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_ova_CNA_sort_common_gene_11859.txt", row.names = 1)
ova_dna_num = ova_dna
indx <- sapply(ova_dna, is.factor)
ova_dna_num[indx] = lapply(ova_dna[indx], function(x) as.numeric(as.character(x)))

ova_jhu_proteome = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_ova_JHU_proteome_sort_common_gene_7061.txt", row.names = 1, fill = TRUE)
ova_jhu_proteome_num = ova_jhu_proteome
indx <- sapply(ova_jhu_proteome, is.factor)
ova_jhu_proteome_num[indx] = lapply(ova_jhu_proteome[indx], function(x) as.numeric(as.character(x)))

ova_pnnl_proteome = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt", row.names = 1, fill = TRUE)
ova_pnnl_proteome_num = ova_pnnl_proteome
indx <- sapply(ova_pnnl_proteome, is.factor)
ova_pnnl_proteome_num[indx] = lapply(ova_pnnl_proteome[indx], function(x) as.numeric(as.character(x)))

ova_rna = read.table("/fs/project/PAS0272/Tara/proteomics_dream_challenge/retrospective_ova_rna_seq_sort_common_gene_15121.txt", row.names = 1, fill = TRUE)
ova_rna_num = ova_rna
indx <- sapply(ova_rna, is.factor)
ova_rna_num[indx] = lapply(ova_rna[indx], function(x) as.numeric(as.character(x)))

#Show box plots for each data set.
png('breast_dna_boxplot.png')
boxplot(t(breast_dna_num)[,2:51])
dev.off()

png('breast_proteome_boxplot.png')
boxplot(t(breast_proteome_num)[,2:51])
dev.off()

png('breast_rna_boxplot.png')
boxplot(t(breast_rna_num)[,2:51])
dev.off()

png('ova_microarray_boxplot.png')
boxplot(t(ova_microarray_num)[,2:51])
dev.off()

png('ova_dna_boxplot.png')
boxplot(t(ova_dna_num)[,2:51])
dev.off()

png('ova_jhu_proteome_boxplot.png')
boxplot(t(ova_jhu_proteome_num)[,2:51])
dev.off()

png('ova_pnnl_proteome_boxplot.png')
boxplot(t(ova_pnnl_proteome_num)[,2:51])
dev.off()

png('ova_rna_boxplot.png')
boxplot(t(ova_rna_num)[,2:51])
dev.off()

#Print summaries for each data set.
sink('breast_dna_summary.txt', append = FALSE)
summary(t(breast_dna))
sink()

sink('breast_proteome_summary.txt', append = FALSE)
summary(t(breast_proteome))
sink()

sink('breast_rna_summary.txt', append = FALSE)
summary(t(breast_rna))
sink()

sink('ova_dna_summary.txt', append = FALSE)
summary(t(ova_dna))
sink()

sink('ova_microarray_summary.txt', append = FALSE)
summary(t(ova_microarray))
sink()

sink('ova_jhu_proteome_summary.txt', append = FALSE)
summary(t(ova_jhu_proteome))
sink()

sink('ova_pnnl_proteome_summary.txt', append = FALSE)
summary(t(ova_pnnl_proteome))
sink()

sink('ova_rna_summary.txt', append = FALSE)
summary(t(ova_rna))
sink()

#Get correlation matrices.
breast_list = list(breast_dna_num, breast_proteome_num, breast_rna_num)
common_names_breast = Reduce(intersect, lapply(breast_list, row.names))
common_patients_breast = Reduce(intersect, lapply(breast_list, colnames))
breast_list = lapply(breast_list, function(x) { x[row.names(x) %in% common_names_breast,] })
breast_list = lapply(breast_list, function(x) { x[,colnames(x) %in% common_patients_breast] })
breast_dna_corr = cor(t(breast_list[[1]][-1,-1]), t(breast_list[[2]][-1,-1]), use = "pairwise.complete.obs")
saveRDS(breast_dna_corr, file = "breast_dna_corr.rds")
breast_rna_corr = cor(t(breast_list[[2]][-1,-1]), t(breast_list[[3]][-1,-1]), use = "pairwise.complete.obs")
saveRDS(breast_rna_corr, file = "breast_rna_corr.rds")

ova_list = list(ova_dna_num, ova_microarray_num, ova_jhu_proteome_num, ova_pnnl_proteome_num, ova_rna_num)
common_names_ova = Reduce(intersect, lapply(ova_list, row.names))
common_patients_ova = Reduce(intersect, lapply(ova_list, colnames))
ova_list = lapply(ova_list, function(x) { x[row.names(x) %in% common_names_ova,] })
ova_list = lapply(ova_list, function(x) { x[,colnames(x) %in% common_patients_ova] })
ova_jhu_dna_corr = cor(t(ova_list[[1]]), t(ova_list[[3]]), use = "pairwise.complete.obs")
saveRDS(ova_jhu_dna_corr, file = "ova_jha_dna_corr.rds")

ova_pnnl_dna_corr = cor(t(ova_list[[1]]), t(ova_list[[4]]), use = "pairwise.complete.obs")
saveRDS(ova_pnnl_dna_corr, file = "ova_pnnl_dna_corr.rds")

ova_jhu_rna_corr = cor(t(ova_list[[5]]), t(ova_list[[3]]), use = "pairwise.complete.obs")
saveRDS(ova_jhu_rna_corr, file = "ova_jhu_rna_corr.rds")

ova_pnnl_rna_corr = cor(t(ova_list[[5]]), t(ova_list[[4]]), use = "pairwise.complete.obs")
saveRDS(ova_pnnl_rna_corr, file = "ova_pnnl_rna_corr.rds")

ova_jhu_microarray_corr = cor(t(ova_list[[2]]), t(ova_list[[3]]), use = "pairwise.complete.obs")
saveRDS(ova_jhu_microarray_corr, file = "ova_jhu_microarray_corr.rds")

ova_pnnl_microarray_corr = cor(t(ova_list[[2]]), t(ova_list[[4]]), use = "pairwise.complete.obs")
saveRDS(ova_pnnl_microarray_corr, file = "ova_pnnl_microarray_corr.rds")


#Show correlation plots.
library(reshape2)
library(ggplot2)

melted_cormat <- melt(breast_dna_corr)
png('breast_dna_corrplot.png')
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
dev.off()