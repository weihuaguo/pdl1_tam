# This script is designed to analyze in-house scRNA-seq data 
# Weihua Guo, Ph.D.
# 02/14/2020
# Cleaned on 09/27/2023

srQCPrePro <- function(srsc, plotPf, res = 300) {
	cat("Start to QC and pre-processing Seurat object...\n")
	mito.genes <- grep(pattern = "^MT-", x = rownames(srsc@assays[["RNA"]]), value = TRUE)
	srsc[["percent.mt"]] <- PercentageFeatureSet(srsc, pattern = "^MT-")

	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6)

	qcsc1 <- FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "percent.mt")
	qcsc2 <- FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	cqcgg <- CombinePlots(plots = list(qcsc1, qcsc2))
	ggsave(plot = cqcgg, filename = paste(plotPf, "QCPlot2.png", sep = ""), dpi = res, width = 9, height = 6)
	srsc <- subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

	###################
	srsc <- NormalizeData(srsc, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

	################
	srsc <- FindVariableFeatures(srsc, selection.method = "vst", nfeatures = 2000)
	top10Features <- head(VariableFeatures(srsc), 18)
	varGG <- VariableFeaturePlot(srsc)
	labVargg <- LabelPoints(plot = varGG, points = top10Features, repel = TRUE)
	ggsave(plot = labVargg, filename = paste(plotPf, "variable_features.png", sep = ""), dpi = res, width = 9, height = 6)
	return(srsc)
}

inteSrPro <- function(srsc, plotPf, batchName = NULL, res = 300, ndim = 18, rdsSave = FALSE, plotFeat = NULL, jsFlag = FALSE, 
		      clusterRes = 0.5) {
	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = batchName)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6) # FigS1x

	srsc <- ScaleData(srsc)

	#############
	srsc <- RunPCA(srsc)
	if (FALSE) { #FPCA
	pcaGG1 = VizDimLoadings(srsc, dim = 1:6, reduction = "pca")
	ggsave(plot = pcaGG1, filename = paste(plotPf, "PCA1.png", sep = ""), dpi = res, width = 9, height = 9)
	
	pcaGG2 = DimPlot(srsc, reduction = "pca")
	ggsave(plot = pcaGG2, filename = paste(plotPf, "PCA2.png", sep = ""), dpi = res, width = 9, height = 6)

	png(paste(plotPf, "PCA_HEATMAP.png", sep = ""), res = res, width = 9, height = 16, units = "in")
	DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
	gar = dev.off()
	} #FPCA

	if (jsFlag) {
		jst = Sys.time()
		cat("Running Jack Straw method to determine the PC dimension to use...\n")
		srsc <- JackStraw(srsc, num.replicate = 100)
		srsc <- ScoreJackStraw(srsc, dims = 1:20)
		jsGG <- JackStrawPlot(srsc, dims = 1:20)
		ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), dpi = res, width = 9, heigh = 6) #FigS1x
		cat("Jack Straw costs")
		print(Sys.time()-jst)
	}

	elGG <- ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6) #FigS1x

	srsc <- FindNeighbors(srsc, dims = 1:ndim)
	srsc <- FindClusters(srsc, resolution = clusterRes)

	srsc <- RunUMAP(srsc, dims = 1:ndim)
	
	umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
	ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)

	if (!is.null(plotFeat)) {
		cat("Start to plot feature plots and violin plots\n")
		vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "counts", log = TRUE, ncol = 3, pt.size = 0.5)
		ggsave(plot = vlnGG, filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30) # FigS1x

		featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
		ggsave(plot = featGG, filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30) #FigS1x

		dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()
		ggsave(plot = dotGG, filename = paste(plotPf, "DotPlot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 9, height = 6) #FigS1x
	}

	if (!is.null(batchName)) {
		batchCheckGG <- DimPlot(srsc, reduction = "umap", group.by = batchName)
		ggsave(plot = batchCheckGG, filename = paste(plotPf, "batch_effect_check_dimPlot.png", sep = ""),
		       dpi = res, width = 9, height = 6) #FigS1x
	}
	
	srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

	write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_MARKERS.csv", sep = "")) #ST1x
	write.csv(topMarkers, file = paste(plotPf, "top10_pos_markers.csv", sep = "")) 

	topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
	png(paste(plotPf, "top5_marker_heatmap.png", sep = ""), res = res, width = 16, height = 18, units = "in")
	print(DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend())
	gar = dev.off()

	topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
	geneMarkers <- unique(as.character(topDotMarkers$gene))
#	geneMarkers = unique(c(as.character(topDotMarkers$gene), "CD14", "CD68", "HLA-DRA", "CD163", "STAB1"))

	dotMarkerGG <- DotPlot(srsc, features = geneMarkers) + RotatedAxis() + coord_flip()
	ggsave(plot = dotMarkerGG, filename = paste(plotPf, "Top5Marker_DotPlot.png", sep = ""), 
	       dpi = res, width = 9, height = 16)

	
	if (rdsSave) {
		saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	}
	print(srsc)
	return(srsc)
}


tst <- Sys.time()
suppressMessages(library(Seurat))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
#suppressMessages(library(MAST))

workDir <- "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub"
dataDir <- paste(workDir, "1_inhouse_scrna_seq", sep = "/") # Please change this to repeat the work
resDir <- paste(workDir, "2_inhouse_scrna_seq_result", sep = "/")

expID <- "inhouse_3prime_tumor_v323"
readFlag <- TRUE
inteFlag <- TRUE
inteFilter <- 100 # Minimum cell number for integration
clstFlag <- TRUE

dir.create(file.path(resDir, expID), showWarnings = FALSE)
expDir <- paste(resDir, expID, sep = "/")

intePf <- paste(expDir, "/", expID, "_Integrated_Results_", sep = "")


if (readFlag) {
	inFiles <- list.files(dataDir, pattern = "_filtered_gene_bc_matrices")
	saName <- str_split_fixed(inFiles, "_fil", n = 2)[,1]

	cat("All the samples:", saName, "\n")
	## Initital a dataframe to record the basic process results
	recColName <- c("tissue", "patient", "before_cell", "after_cell", "before_nFeature", "after_nFeature", 
			"median_mtper", "mean_mtper", "cell_kept")
	inteRecDf <- as.data.frame(matrix(ncol = length(recColName), nrow = length(saName)))
	rownames(inteRecDf) <- saName
	colnames(inteRecDf) <- recColName
	## Initial list for integration
	srData <- vector(mode = "list", length = length(saName))
	srRawObjs <- vector(mode = "list", length = length(saName))
	srCleanObjs <- vector(mode = "list", length = length(saName))
	names(srData) <- saName
	names(srRawObjs) <- saName
	names(srCleanObjs) <- saName

	for (isa in 1:length(inFiles)) {
		st <- Sys.time()
		cat(saName[isa], ":" , inFiles[isa], "\n")
		tmpDir <- paste(dataDir, inFiles[isa], sep = "/")
		tmpData <- try(Read10X(tmpDir))
		if (class(tmpData) == "try-error") {
			stop("STILL NOT WORKING!!!")
		}
		colnames(tmpData) <- str_c(saName[isa], "_", colnames(tmpData))
#		print(head(colnames(tmpData)))
		tmpData <- as(tmpData, 'CsparseMatrix')
#		print(class(tmpData))

		## NOTE: need some basic filter to remove totally empty droplet!!!
		tmpObj <- CreateSeuratObject(counts = tmpData, project = saName[isa], min.cells = 3, min.features = 5)
		print(tmpObj)
#		q(save = "no")

		tisName <- str_split_fixed(saName[isa], "_", n = 3)[,3]
		rnaDirc <- str_split_fixed(saName[isa], "_", n = 3)[,2]
		patName <- str_split_fixed(saName[isa], "_", n = 3)[,1]
		cat(tisName, patName, rnaDirc, "\n")
		tmpObj@meta.data$tissue <- tisName
		tmpObj@meta.data$patient <- patName
		tmpObj@meta.data$sampleID <- saName[isa]

		srData[saName[isa]] <- tmpData
		srRawObjs[saName[isa]] <- tmpObj

		plotPrefix <- paste(expDir, "/", expID, "_", saName[isa], "_", sep = "")
		tmpProObj <- srQCPrePro(tmpObj, plotPf = plotPrefix)
		srCleanObjs[saName[isa]] <- tmpProObj

		inteRecDf[saName[isa], "tissue"] <- tisName
		inteRecDf[saName[isa], "patient"] <- patName
		inteRecDf[saName[isa], "before_cell"] <- length(colnames(tmpObj))
		inteRecDf[saName[isa], "after_cell"] <- length(colnames(tmpProObj))
		inteRecDf[saName[isa], "before_nFeature"] <- median(tmpObj@meta.data$nFeature_RNA)
		inteRecDf[saName[isa], "after_nFeature"] <- median(tmpProObj@meta.data$nFeature_RNA)
		inteRecDf[saName[isa], "median_mtper"] <- median(tmpProObj@meta.data$percent.mt)
		inteRecDf[saName[isa], "mean_mtper"] <- mean(tmpProObj@meta.data$percent.mt)

		print(Sys.time() - st)
		cat("+++++++++++++++++++++++++++++++++++++++++++++\n\n")
	}
	inteRecDf$cell_kept <- inteRecDf$after_cell/inteRecDf$before_cell
	write.csv(inteRecDf, paste(expDir, "/", expID, "_integration_records.csv", sep = ""))

	cat("Saving Seurat object list to a RDS file...\n")
	dataRDS <- paste(expDir, "/", expID, "_all_filtered_10X_data_list.RDS", sep = "")
	rawObjRDS <- paste(expDir, "/", expID, "_all_raw_filtered_10X_seurat_objects_list.RDS", sep = "")
	cleanObjRDS <- paste(expDir, "/", expID, "_all_clean_filtered_10X_seurat_objects_list.RDS", sep = "")
	saveRDS(srData, file = dataRDS)
	saveRDS(srRawObjs, file = rawObjRDS)
	saveRDS(srCleanObjs, file = cleanObjRDS)
}

if (inteFlag) {
	srCleanObjs <- readRDS(paste(expDir, "/", expID, "_all_clean_filtered_10X_seurat_objects_list.RDS", sep = ""))
	st <- Sys.time()
	cat("Start to merge the samples...\n")
	integrateList <- srCleanObjs
	for (inm in names(srCleanObjs)) {
		if(length(colnames(srCleanObjs[[inm]])) <= inteFilter) {integrateList[[inm]] <- NULL}
	}
	integrateAnchors <- FindIntegrationAnchors(object.list = integrateList, dims = 1:50, k.filter = inteFilter)
	inteObjs <- IntegrateData(anchorset = integrateAnchors, dims = 1:50)
	DefaultAssay(inteObjs) <- "integrated"
	print(inteObjs)
	print(Sys.time() - st)

	saveRDS(inteObjs, file = paste(intePf, "unprocessed_integrated_seurat_object.RDS", sep = ""))
}

if (clstFlag) {
	cat("Analyzing integrated Seurat object...\n")
	sigMarkers <- c("EPCAM", "KRT19", "COL1A1", "DCN", "VWF", "PECAM1", "CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA-DRA", "NCAM1", "KIT", "FCER1A")
	pst <- Sys.time()
	inteObjs <- readRDS(paste(intePf, "unprocessed_integrated_seurat_object.RDS", sep = ""))

	inteProObj <- inteSrPro(inteObjs, batchName = c("patient"), plotPf = intePf, rdsSave = TRUE, ndim = 20, plotFeat = sigMarkers, clusterRes=0.5)
	cat("Total")
	print(Sys.time()-pst)
}

# TODO: Annotated cell types (UMAP, violin of key markers, heatmap, dotplot, cell proportion crossing patients, cancer TAM ratio, batch effect)
# TODO: Subcluster TAM/TAM+DC, check the PD-L1+ and SIGLEC15+ cell number, screen the resolution
# TODO: Build a clustree to determine that optimal resolution
# TODO: TAM subcluster: UMAP, marker heatmap and plots
# TODO: Supervised clustering????
# TODO: CCI
