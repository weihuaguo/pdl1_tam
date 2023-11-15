rm(list = ls())

tst <- Sys.time()
suppressMessages(library(Seurat))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(MAST))
suppressMessages(library(clustree))

outDir <- "/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020"
resObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res020_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res020_Seurat_Objects_Clustered.RDS")
inteProObj <- resObj
colnames(inteProObj@meta.data)[colnames(inteProObj@meta.data) == "integrated_snn_res.0.5"] <- "integrated_snn_res.0.0"
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res030_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res030_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.3 <- nxtObj@meta.data$integrated_snn_res.0.3
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res040_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res040_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.4 <- nxtObj@meta.data$integrated_snn_res.0.4
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res050_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res050_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.5 <- nxtObj@meta.data$integrated_snn_res.0.5
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res080_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res080_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.8 <- nxtObj@meta.data$integrated_snn_res.0.8
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res100_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res100_Seurat_Objects_Clustered.RDS")

inteProObj@meta.data$integrated_snn_res.1.0 <- nxtObj@meta.data$integrated_snn_res.1

cltr <- clustree(inteProObj, prefix = "integrated_snn_res.") + scale_color_brewer(palette = "Set1")
ggsave(paste(outDir, "/clean3prime_tumor_02262020_clustree.png", sep = ""), cltr, dpi = 300, width = 9, height = 9)


print(head(inteProObj@meta.data))
q(save = "no")


allObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub/2_inhouse_scrna_seq_result/inhouse_3prime_tumor_final/myeloid_sub_cluster/tam_subcluster_res_1.2/tam_subcluster_res_1.2_Seurat_Objects_Clustered.RDS")
oldObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res050_subtype_analysis/clean3prime_tumor_02262020_Macrophage_Res050_Seurat_Objects_Clustered.RDS")
outDir <- "/home/weihua/mnts/smb_plee/Group/weihua/sc_old_result/AllThreePrimeRNAResults/clean3prime_tumor_02262020/Macrophage_Res050_subtype_analysis"

print(head(allObj@meta.data))
print(head(oldObj@meta.data))
ress <- c(seq(0.1,1.2, 0.1))#, seq(1.0,2.0,0.2))
allMeta <- merge(allObj@meta.data, oldObj@meta.data, by = c("row.names", "orig.ident"), suffixes = c("", "_old"))
print(head(allMeta))
ctsNewOld <- allMeta %>%
	group_by(`integrated_snn_res.0.5`, `integrated_snn_res.0.5_old`) %>%
	summarize(n = n())
ctsNewOld <- ctsNewOld %>%
	group_by(`integrated_snn_res.0.5_old`) %>%
	mutate(nOld = sum(n))
ctsNewOld$rOld <- ctsNewOld$n/ctsNewOld$nOld*100
tile_gg <- ggplot(ctsNewOld, aes(x = `integrated_snn_res.0.5_old`, y = `integrated_snn_res.0.5`, fill = rOld, 
				 label = sprintf("%0.2f", round(rOld, digits = 2)))) +
	geom_tile()+
	geom_text()+
	scale_fill_distiller(palette = "Spectral") +
	theme_classic()
ggsave(paste(outDir, "/old_new_res0.5_heatmap.png", sep = ""), dpi = 300, width = 6, height = 6)

inteProObj <- allObj
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 3] <- 4
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 4] <- 3
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 5] <- 3
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 6] <- 5
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 7] <- 6
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 8] <- 7
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 9] <- 8
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 11] <- 9
inteProObj@meta.data$integrated_snn_res.0.5[inteProObj@meta.data$integrated_snn_res.0.5 == 12] <- 11

saveRDS(inteProObj, "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub/2_inhouse_scrna_seq_result/inhouse_3prime_tumor_final/myeloid_sub_cluster/tam_subcluster_res_1.2/tam_subcluster_res_1.2_Seurat_Objects_Corrected_Clustered.RDS")
all_cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(ress)+1)
cltr <- clustree(inteProObj, prefix = "integrated_snn_res.") + scale_color_manual(values = all_cols)
ggsave(paste(outDir, "/corrected_clustree.png", sep = ""), cltr, dpi = 300, width = 12, height = 12)

