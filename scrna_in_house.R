# This script is designed to analyze in-house scRNA-seq data 
# Weihua Guo, Ph.D.
# 02/14/2020
# Cleaned on 09/27/2023

rm(list = ls())
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

inteSrPro <- function(srsc, plotPf, batchName = NULL, res = 300, ndim = 18, rdsSave = FALSE, plotFeat = NULL, jsFlag = FALSE) {
	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = batchName)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6) 

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
		ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), dpi = res, width = 9, heigh = 6)
		cat("Jack Straw costs")
		print(Sys.time()-jst)
	}

	elGG <- ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)

	srsc <- FindNeighbors(srsc, dims = 1:ndim)
	srsc <- FindClusters(srsc)

	srsc <- RunUMAP(srsc, dims = 1:ndim)
	
	umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
	ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)

	if (!is.null(plotFeat)) {
		cat("Start to plot feature plots and violin plots\n")
		vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "data", log = TRUE, ncol = 3, pt.size = 0.5)
		ggsave(plot = vlnGG, filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30) 

		featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
		ggsave(plot = featGG, filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30)

		dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()
		ggsave(plot = dotGG, filename = paste(plotPf, "DotPlot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}

	if (!is.null(batchName)) {
		batchCheckGG <- DimPlot(srsc, reduction = "umap", group.by = batchName)
		ggsave(plot = batchCheckGG, filename = paste(plotPf, "batch_effect_check_dimPlot.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}
	
	srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

	write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_MARKERS.csv", sep = "")) #ST1x
	write.csv(topMarkers, file = paste(plotPf, "top10_pos_markers.csv", sep = "")) 

	topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	hm <- DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend()
	ggsave(paste(plotPf, "top5_marker_heatmap.png", sep = ""), hm, dpi = res, width = 24, height = 18)

	topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	geneMarkers <- unique(as.character(topDotMarkers$gene))

	dotMarkerGG <- DotPlot(srsc, features = geneMarkers) + RotatedAxis() + coord_flip()
	ggsave(plot = dotMarkerGG, filename = paste(plotPf, "Top5Marker_DotPlot.png", sep = ""), 
	       dpi = res, width = 9, height = 16)

	
	if (rdsSave) {
		saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	}
	print(srsc)
	return(srsc)
}

resScreen <- function(srsc, plotPf, plotDir, batchName = NULL, res = 300, ndim = 18, rdsSave = FALSE, plotFeat = NULL, jsFlag = FALSE, resolutions) {
	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = batchName)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6) 

	srsc <- ScaleData(srsc)

	#############
	srsc <- RunPCA(srsc)

	elGG <- ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)
	srsc <- FindNeighbors(srsc, dims = 1:ndim)
	srsc <- RunUMAP(srsc, dims = 1:ndim)

	if (!is.null(batchName)) {
		batchCheckGG <- DimPlot(srsc, reduction = "umap", group.by = batchName)
		ggsave(plot = batchCheckGG, filename = paste(plotPf, "batch_effect_check_dimPlot.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}


	for (ires in resolutions) {
		resFolder <- paste("tam_subcluster_res", ires, sep = "_")
		dir.create(file.path(plotDir, resFolder), showWarnings = FALSE)
		resDir <- paste(plotDir, resFolder, sep = "/")
		rppf <- paste(resDir, "/", resFolder, "_", sep = "")

		srsc <- FindClusters(srsc, resolution = ires)

		srsc <- BuildClusterTree(srsc)
		png(paste(rppf, "cluster_tree.png", sep = ""), res = 300, width = 12, height = 18, unit = 'in')
		PlotClusterTree(object = srsc)
		gar <- dev.off()

		umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
		ggsave(plot = umapGG, filename = paste(rppf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)

		if (!is.null(plotFeat)) {
			cat("Start to plot feature plots and violin plots\n")
			vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "counts", log = TRUE, ncol = 3, pt.size = 0.5)
			ggsave(plot = vlnGG, filename = paste(rppf, "violin_plot_for_cell_annotation.png", sep = ""), 
			       dpi = res, width = 16, height = 16) 

			featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
			ggsave(plot = featGG, filename = paste(rppf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
			       dpi = res, width = 16, height = 16)

			dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()
			ggsave(plot = dotGG, filename = paste(rppf, "DotPlot_for_cell_annotation.png", sep = ""),
			       dpi = res, width = 9, height = 6)
		}

		
		srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

		write.csv(srsc.markers, file = paste(rppf, "ALL_POS_MARKERS.csv", sep = "")) #ST1x
		write.csv(topMarkers, file = paste(rppf, "top10_pos_markers.csv", sep = "")) 

		topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
		hm <- DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend()
		ggsave(paste(rppf, "top5_marker_heatmap.png", sep = ""), hm, dpi = res, width = 24, height = 18)

		topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
		geneMarkers <- unique(as.character(topDotMarkers$gene))

		dotMarkerGG <- DotPlot(srsc, features = geneMarkers) + RotatedAxis() + coord_flip()
		ggsave(plot = dotMarkerGG, filename = paste(rppf, "Top5Marker_DotPlot.png", sep = ""), 
		       dpi = res, width = 9, height = 16)

		
		if (rdsSave) {
			saveRDS(srsc, file = paste(rppf, "Seurat_Objects_Clustered.RDS", sep = ""))
		}
		print(srsc)
	}
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
suppressMessages(library(RColorBrewer))
suppressMessages(library(MAST))
suppressMessages(library(clustree))


workDir <- "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub"
dataDir <- paste(workDir, "1_inhouse_scrna_seq", sep = "/") # Please change this to repeat the work
resDir <- paste(workDir, "2_inhouse_scrna_seq_result", sep = "/")

expID <- "inhouse_3prime_tumor"
readFlag <- FALSE
inteFlag <- FALSE
inteFilter <- 100 # Minimum cell number for integration
clstFlag <- FALSE
majorVisFlag <- FALSE
tamClusterFlag <- TRUE

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

		## NOTE: need some basic filter to remove totally empty droplet!!!
		tmpObj <- CreateSeuratObject(counts = tmpData, project = saName[isa], min.cells = 3, min.features = 5)
		print(tmpObj)

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

	inteProObj <- inteSrPro(inteObjs, batchName = c("patient"), plotPf = intePf, rdsSave = TRUE, ndim = 20, plotFeat = sigMarkers)
	cat("Total")
	print(Sys.time()-pst)
}

if (majorVisFlag) {
	cat("Merge and visualize all the cells\n")
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))
	cellAnnDf <- as.data.frame(read_excel(paste(intePf, "ALL_POS_MARKERS.xlsx", sep = "")))
	rownames(cellAnnDf) <- cellAnnDf[,1]
	cellAnnDf[,1] <- NULL
	cellMarkerDf <- cellAnnDf
	cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$cell_type),]
	print(head(cellAnnDf))
	proObj@meta.data$major_cell_type <- "NA"
	for (ic in cellAnnDf$cluster) {
		proObj@meta.data$major_cell_type[proObj@meta.data$seurat_clusters == ic] <- cellAnnDf$cell_type[cellAnnDf$cluster == ic]
	}
	print(head(proObj@meta.data))

	sigMarkers <- c("EPCAM", "KRT19", "COL1A1", "DCN", "VWF", "PECAM1", "CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA-DRA", "CD1C", "KIT", "ITGAM")
	sigExpr <- proObj[['integrated']]@data[sigMarkers,]
	sigExpr <- as.data.frame(sigExpr)
	sigExpr$gene <- rownames(sigExpr)
#	print(dim(sigExpr))
#	print(sigExpr[1:9,1:6])

	gathExpr <- gather(sigExpr, "cell_id", "expr", rownames(proObj@meta.data))
	gathExpr <- merge(gathExpr, proObj@meta.data, by.x = "cell_id", by.y = "row.names", all.x = TRUE)
	gathExpr$gene <- factor(gathExpr$gene, levels = sigMarkers)
#	print(head(gathExpr))

	figs1e <- ggplot(gathExpr, aes(x = seurat_clusters, y = expr, fill = seurat_clusters)) +
		geom_violin(scale = "width", trim = TRUE) +
		facet_wrap(.~gene, nrow = 1, scales = "free_x") +
		labs(x = "Cluster (r = 0.8)", y = "Normalized expression", fill = "Cluster (r = 0.8)") +
		coord_flip() +
		theme_classic()
	ggsave(plot = figs1e, filename = paste(intePf, "FigS1E_sc_violin_plot_for_cell_annotation.png", sep = ""), 
	       dpi = 300, width = 16, height = 4.5) # FigS1x
	saveRDS(figs1e, paste(intePf, "FigS1E_sc_violin_plot_for_cell_annotation.RDS", sep = ""))


	figs1j <- ggplot(gathExpr, aes(x = major_cell_type, y = expr, fill = major_cell_type)) +
		geom_violin(scale = "width", trim = TRUE) +
		facet_wrap(.~gene, nrow = 1, scales = "free_x") +
		scale_fill_brewer(palette = "Spectral") +
		labs(x = "General cell type in TME", y = "Normalized expression", fill = "General cell type\nin TME") +
		coord_flip() +
		theme_classic()
	ggsave(plot = figs1j, filename = paste(intePf, "FigS1J_mc_violin_plot_for_cell_annotation.png", sep = ""), 
	       dpi = 300, width = 16, height = 2.5) # FigS1x
	saveRDS(figs1j, paste(intePf, "FigS1J_mc_violin_plot_for_cell_annotation.RDS", sep = ""))

	cat("Heatmaps\n")
	topMarkers <- cellMarkerDf %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	figs1d <- DoHeatmap(object = proObj, features = topMarkers$gene) + NoLegend()
	ggsave(paste(intePf, "FigS1D_top5_marker_heatmap.png", sep = ""), figs1d, dpi = 300, width = 24, height = 18)
	saveRDS(figs1d, paste(intePf, "FigS1D_sc_top5_marker_heatmap.RDS", sep = ""))

	if (FALSE) { # NOTE: To save some time
	c <- 0
	for (imct in unique(proObj@meta.data$major_cell_type)) {
		cat("Caculating", imct, "\n")
		tmpMarkers <- FindMarkers(proObj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, group.by = "major_cell_type", ident.1 = imct)
		tmpMarkers$gene <- rownames(tmpMarkers)
		tmpMarkers$cluster <- imct
		if (c == 0) {
			gnrlMarkerDf <- tmpMarkers
		} else {
			gnrlMarkerDf <- rbind(gnrlMarkerDf, tmpMarkers)
		}
		c <- c+1
	}
	print(head(gnrlMarkerDf))
	write.csv(gnrlMarkerDf, paste(intePf, "ALL_POS_MARKERS_MAJOR_CELL_TYPE.csv"))
	} else {
		gnrlMarkerDf <- read.csv(paste(intePf, "ALL_POS_MARKERS_MAJOR_CELL_TYPE.csv"), row.names = 1, header = T, check.names = F)
	}

	cols <- brewer.pal(length(unique(gnrlMarkerDf$cluster)), "Spectral")

	topMarkers <- gnrlMarkerDf %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
	topOrders <- topMarkers$gene[order(topMarkers$cluster, topMarkers$avg_log2FC, decreasing = T)]
	topOrders <- topOrders[!duplicated(topOrders)]
	figs1i <- DoHeatmap(object = proObj, features = topMarkers$gene, group.colors = cols, group.by = "major_cell_type") + NoLegend()
#	print(figs1i$data)
	figs1i$data$Feature <- factor(figs1i$data$Feature, levels = topOrders)
	ggsave(paste(intePf, "FigS1I_top10_marker_mct_heatmap.png", sep = ""), figs1i, dpi = 300, width = 24, height = 18)
	saveRDS(figs1i, paste(intePf, "FigS1I_mc_top10_marker_heatmap.RDS", sep = ""))

	cat("Cell proportions\n")
	scCts <- proObj@meta.data %>%
		group_by(patient, seurat_clusters) %>%
		summarize(n = n())

	figs1c <- ggplot(scCts, aes(x = patient, y = n, fill = seurat_clusters)) +
		geom_bar(stat = "identity", position = "fill", color = "gray") +
		coord_flip() +
		labs(x = "Sample ID", y = "Cell proportion", fill = "Cluster (r = 0.8)") +
		theme_classic()
	ggsave(plot = figs1c, filename = paste(intePf, "FigS1C_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1c, paste(intePf, "FigS1C.RDS", sep = ""))

	mcCts <- proObj@meta.data %>%
		group_by(patient, major_cell_type) %>%
		summarize(n = n())
	print(head(mcCts))

	figs1h <- ggplot(mcCts, aes(x = patient, y = n, fill = major_cell_type)) +
		geom_bar(stat = "identity", position = "fill", color = "gray") +
		scale_fill_brewer(palette = "Spectral") +
		coord_flip() +
		labs(x = "Sample ID", y = "Cell proportion", fill = "General cell type\nin TME") +
		theme_classic()
	ggsave(plot = figs1h, filename = paste(intePf, "FigS1H_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1h, paste(intePf, "FigS1H.RDS", sep = ""))

	proObj@meta.data$PDL1_expr <- proObj[["RNA"]]@data["CD274",]
	proObj@meta.data$PDL1_status <- ifelse(proObj@meta.data$PDL1_expr > 0, "PDL1+", "PDL1-")

	pdl1Cts <- proObj@meta.data %>%
		filter(major_cell_type == "Cancer cell" | major_cell_type == "TAM") %>%
		group_by(patient, major_cell_type, PDL1_status) %>%
		summarize(n = n())

	pdl1Spr <- as.data.frame(spread(pdl1Cts, "major_cell_type", "n"))
	pdl1Spr$CM_ratio <- pdl1Spr$`Cancer cell`/pdl1Spr$TAM

	cmSpr <- mcCts %>%
		filter(major_cell_type == "Cancer cell"|major_cell_type == "TAM") %>%
		spread(major_cell_type, n)
	cmSpr$CM_ratio <- cmSpr$`Cancer cell`/cmSpr$TAM
	cmSpr$PDL1_status <- "PDL1+ & PDL1-"

	cmRatio <- rbind(cmSpr, pdl1Spr)

	figs1g <- ggplot(cmRatio, aes(x = PDL1_status, y = CM_ratio, fill = patient)) +
		geom_bar(stat = "identity", color = "gray", position = "dodge") +
		scale_fill_brewer(palette = "Dark2") +
		coord_flip() +
		labs(x = "PD-L1 status", y = "Cancer cell/TAM ratio", fill = "Sample ID") +
		theme_classic()
	ggsave(plot = figs1g, filename = paste(intePf, "FigS1G_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1g, paste(intePf, "FigS1G.RDS", sep = ""))

	cat("UMAPs\n")
	figs1a <- DimPlot(proObj, reduction = "umap", label = TRUE) + labs(color = "Clusters (r = 0.8)", title = "")
	ggsave(plot = figs1a, filename = paste(intePf, "FigS1A_UMAP.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1a, paste(intePf, "FigS1A.RDS", sep = ""))

	figs1b <- DimPlot(proObj, reduction = "umap", label = FALSE, group.by = "patient") + 
		labs(color = "Sample ID", title = "Batch effect crossing samples") +
		scale_color_brewer(palette = "Dark2")
	ggsave(plot = figs1b, filename = paste(intePf, "FigS1B_UMAP.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1b, paste(intePf, "FigS1B.RDS", sep = ""))

	figs1f <- DimPlot(proObj, reduction = "umap", label = TRUE, group.by = "major_cell_type") + 
		labs(color = "General cell type\nin TME", title = "Annotated and merged cell types") +
		scale_color_brewer(palette = "Spectral")

	ggsave(plot = figs1f, filename = paste(intePf, "FigS1F_UMAP.png", sep = ""), dpi = 600, width = 9, heigh = 6)
	saveRDS(figs1f, paste(intePf, "FigS1F.RDS", sep = ""))

}

if (tamClusterFlag) {
	cat("Cluster TAMs...\n")
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))

	ist <- Sys.time()
	subSrsc <- subset(proObj, idents = c(4,15,16)) # NOTE: ALL_POS_MARKER.xlsx

	newSrsc <- CreateSeuratObject(counts = subSrsc[["RNA"]]@counts, meta.data = subSrsc@meta.data)
	seuratList <- SplitObject(newSrsc, split.by = "patient")
	
	subFolder <- "tam_sub_cluster"
	dir.create(file.path(expDir, subFolder), showWarnings = FALSE)
	subDir <- paste(expDir, subFolder, sep = "/")
	subPrefix <- paste(subDir, "/", expID, "_", subFolder, "_", sep = "")

	for (i in names(seuratList)) {
		cat("Patient", i, ":", length(colnames(seuratList[[i]])), "\n")
		seuratList[[i]] <- NormalizeData(seuratList[[i]], verbose = FALSE)
		seuratList[[i]] <- FindVariableFeatures(seuratList[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
		if (length(colnames(seuratList[[i]])) <= 50) {
			## NOTE: remove the patients with few cells (<100 cells) for integration!
			seuratList[[i]] <- NULL
			cat("\tPatient", i, "is removed due to small cell number!!!\n")
		}
	}
	integrateAnchors <- FindIntegrationAnchors(object.list = seuratList, k.filter = 50)
	inteObjs <- IntegrateData(anchorset = integrateAnchors, k.weight = 50)
	DefaultAssay(inteObjs) <- "integrated"
	saveRDS(inteObjs, paste(subPrefix, "_Integrated_results.RDS", sep = ""))
	cat("Integration time cost:")
	print(Sys.time()-ist)

	pst <- Sys.time()
	ress <- c(seq(0.1,0.9, 0.1), seq(1.0,2.0,0.2))
	geneOi <- c("CD14", "CD68", "HLA-DRA", "CD274", "SIGLEC15", "PDCD1LG2", "CD1C", "ITGAX", "FCER1A", "TPSAB1", "KIT", "FUT4")
	inteProObj <- resScreen(inteObjs, plotPf = subPrefix, plotDir = subDir, rdsSave = TRUE, jsFlag = FALSE, ndim = 20,
				batchName = c("patient"), plotFeat = geneOi, resolutions = ress)
	cat("Analysis time cost:")
	print(Sys.time()-pst)

	all_cols <- brewer.pal(length(ress), 'Set1')
	colnames(inteProObj@meta.data)[colnames(inteProObj@meta.data) == "integrated_snn_res.1"] <- "integrated_snn_res.0.0"
	cltr <- clustree(inteProObj, prefix = "integrated_snn_res.") + scale_color_manual(values = all_cols)
	ggsave(paste(subPrefix, "clustree.png", sep = ""), cltr, dpi = 300, width = 12, height = 36)

	cat("++++++++++++++++++++++++++++++++++++++++++++++++\n\n")

}
# TODO: Subcluster TAM/TAM+DC, check the PD-L1+ and SIGLEC15+ cell number, screen the resolution
# TODO: Build a clustree to determine that optimal resolution
# TODO: TAM subcluster: UMAP, marker heatmap and plots
# TODO: Supervised clustering????
# TODO: CCI
