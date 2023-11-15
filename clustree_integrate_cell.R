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
outDir <- "/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5"
resObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5/r0.2/r0.2_GSE114727_TUMOR_std_ER_111521_tam_glbr0.5_Seurat_Objects_Clustered.RDS")
inteProObj <- resObj
colnames(inteProObj@meta.data)[colnames(inteProObj@meta.data) == "integrated_snn_res.0.5"] <- "integrated_snn_res.0.0"
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5/r0.4/r0.4_GSE114727_TUMOR_std_ER_111521_tam_glbr0.5_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.4 <- nxtObj@meta.data$integrated_snn_res.0.4
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5/r0.6/r0.6_GSE114727_TUMOR_std_ER_111521_tam_glbr0.5_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.6 <- nxtObj@meta.data$integrated_snn_res.0.6
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5/r0.8/r0.8_GSE114727_TUMOR_std_ER_111521_tam_glbr0.5_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.0.8 <- nxtObj@meta.data$integrated_snn_res.0.8
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5/r1/r1_GSE114727_TUMOR_std_ER_111521_tam_glbr0.5_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.1.0 <- nxtObj@meta.data$integrated_snn_res.1
nxtObj <- readRDS("/home/weihua/mnts/smb_plee/Group/weihua/sc_public_data/breast_results/GSE114727_TUMOR_std_ER_111521/tam_glbr0.5/r1.2/r1.2_GSE114727_TUMOR_std_ER_111521_tam_glbr0.5_Seurat_Objects_Clustered.RDS")
inteProObj@meta.data$integrated_snn_res.1.2 <- nxtObj@meta.data$integrated_snn_res.1.2


cltr <- clustree(inteProObj, prefix = "integrated_snn_res.") + scale_color_brewer(palette = "Set1")
ggsave(paste(outDir, "/clean3prime_tumor_02262020_clustree.png", sep = ""), cltr, dpi = 300, width = 9, height = 9)


print(head(inteProObj@meta.data))
q(save = "no")

