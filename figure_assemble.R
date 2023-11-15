# Assemble all the figures in R
# Weihua Guo, Ph.D.
# 23.10.05

suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(scales))
suppressMessages(library(cowplot))

sumDir <- "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub/9_RDS"

figres <- 300
figs1_flag <- FALSE
figs3_flag <- TRUE

if (figs3_flag) {
	all_files <- list.files(sumDir)
	figs3a_file <- all_files[str_detect(all_files, "FigS3A")]
	figs3b_file <- all_files[str_detect(all_files, "FigS3B")]
	figs3c_file <- all_files[str_detect(all_files, "FigS3C")]
	figs3d_file <- all_files[str_detect(all_files, "FigS3D")]
	figs3e_file <- all_files[str_detect(all_files, "FigS3E")]
	figs3f_file <- all_files[str_detect(all_files, "FigS3F")]
#	figs3g_file <- all_files[str_detect(all_files, "FigS3G")]
#	figs3h_file <- all_files[str_detect(all_files, "FigS3H")]
#	figs3i_file <- all_files[str_detect(all_files, "FigS3I")]
#	figs3j_file <- all_files[str_detect(all_files, "FigS3J")]

	figs3a <- readRDS(paste(sumDir, figs3a_file, sep = "/"))
	figs3b <- readRDS(paste(sumDir, figs3b_file, sep = "/"))
	figs3c <- readRDS(paste(sumDir, figs3c_file, sep = "/"))
	figs3d <- readRDS(paste(sumDir, figs3d_file, sep = "/"))
	figs3e <- readRDS(paste(sumDir, figs3e_file, sep = "/"))
	figs3f <- readRDS(paste(sumDir, figs3f_file, sep = "/"))
#	figs3g <- readRDS(paste(sumDir, figs3g_file, sep = "/"))
#	figs3h <- readRDS(paste(sumDir, figs3h_file, sep = "/"))
#	figs3i <- readRDS(paste(sumDir, figs3i_file, sep = "/"))
#	figs3j <- readRDS(paste(sumDir, figs3j_file, sep = "/"))


	figs3a <- figs3a + theme(text = element_text(size = 16), legend.position = "none")
	figs3b <- figs3b + theme(text = element_text(size = 16))
	figs3c <- figs3c + theme(text = element_text(size = 16), legend.position = "none")
	figs3d <- figs3d + theme(text = element_text(size = 16), axis.text.x = element_blank())
	figs3e <- figs3e + theme(text = element_text(size = 16), axis.text.x = element_blank())
	figs3f <- figs3f + theme(text = element_text(size = 16))
#	figs3g <- figs3g + theme(text = element_text(size = 16))
#	figs3h <- figs3h + theme(text = element_text(size = 16))
#	figs3i <- figs3i + theme(text = element_text(size = 16))
#	figs3j <- figs3j + theme(text = element_text(size = 16))

	upper <- ggarrange(
			   ggarrange(
				     ggarrange(figs3a, figs3b, labels = c("A", "B"), ncol = 2, nrow = 1, font.label = list(size = 20)), 
				     figs3d, figs3e, labels = c(NA, "D", "E"), ncol = 1, nrow = 3, heights = c(1,1.5,1), font.label = list(size = 20), common.legend = F),
			   figs3c, labels = c(NA, "C"), ncol = 2, nrow = 1, widths = c(3,2.4), font.label = list(size = 20))

#	lower <- ggarrange(figs3f, ggarrange(figs3g, figs3h, labels = c("G", "H"), ncol = 2, nrow = 1, font.label = list(size = 20)), 
#			   labels = c("F", NA), ncol = 1, nrow = 2, heights = c(1,1.5), font.label = list(size = 20))

	test <- ggarrange(upper, figs3f, nrow = 2, ncol = 1, labels = c(NA, "F"), heights = c(3,1))
	png(paste(sumDir, "figs3_cowplot.png", sep = "/"), res = figres, width = 32, height = 18, units = 'in')
	print(test)
	gar <- dev.off()

}


if (figs1_flag) {
	all_files <- list.files(sumDir)
	figs1a_file <- all_files[str_detect(all_files, "FigS1A")]
	figs1b_file <- all_files[str_detect(all_files, "FigS1B")]
	figs1c_file <- all_files[str_detect(all_files, "FigS1C")]
	figs1d_file <- all_files[str_detect(all_files, "FigS1D")]
	figs1e_file <- all_files[str_detect(all_files, "FigS1E")]
	figs1f_file <- all_files[str_detect(all_files, "FigS1F")]
	figs1g_file <- all_files[str_detect(all_files, "FigS1G")]
	figs1h_file <- all_files[str_detect(all_files, "FigS1H")]
	figs1i_file <- all_files[str_detect(all_files, "FigS1I")]
#	figs1j_file <- all_files[str_detect(all_files, "FigS1J")]

	figs1a <- readRDS(paste(sumDir, figs1a_file, sep = "/"))
	figs1b <- readRDS(paste(sumDir, figs1b_file, sep = "/"))
	figs1c <- readRDS(paste(sumDir, figs1c_file, sep = "/"))
	figs1d <- readRDS(paste(sumDir, figs1d_file, sep = "/"))
	figs1e <- readRDS(paste(sumDir, figs1e_file, sep = "/"))
	figs1f <- readRDS(paste(sumDir, figs1f_file, sep = "/"))
	figs1g <- readRDS(paste(sumDir, figs1g_file, sep = "/"))
	figs1h <- readRDS(paste(sumDir, figs1h_file, sep = "/"))
	figs1i <- readRDS(paste(sumDir, figs1i_file, sep = "/"))
#	figs1j <- readRDS(paste(sumDir, figs1j_file, sep = "/"))


	figs1a <- figs1a + theme(text = element_text(size = 16))
	figs1b <- figs1b + theme(text = element_text(size = 16))
	figs1c <- figs1c + theme(text = element_text(size = 16), legend.position = "none")
	figs1d <- figs1d + theme(text = element_text(size = 16))
	figs1e <- figs1e + theme(text = element_text(size = 16))
	figs1f <- figs1f + theme(text = element_text(size = 16))
	figs1g <- figs1g + theme(text = element_text(size = 16))
	figs1h <- figs1h + theme(text = element_text(size = 16))
	figs1i <- figs1i + theme(text = element_text(size = 16))
#	figs1j <- figs1j + theme(text = element_text(size = 16))

	upper <- ggarrange(
			   ggarrange(
				     ggarrange(figs1a, figs1b, labels = c("A", "B"), ncol = 2, nrow = 1, font.label = list(size = 20)), 
				     figs1c, figs1e, labels = c(NA, "C", "E"), ncol = 1, nrow = 3, heights = c(1,0.6,2), font.label = list(size = 20), common.legend = F),
			   figs1d, labels = c(NA, "D"), ncol = 2, nrow = 1, widths = c(3,2.4), font.label = list(size = 20))

	lower <- ggarrange(
			   ggarrange(
				     ggarrange(figs1f, figs1g, labels = c("F", "G"), ncol = 2, nrow = 1, font.label = list(size = 20)), 
				     figs1i, labels = c(NA, "I"), ncol = 1, nrow = 2, heights = c(1,1.5), font.label = list(size = 20)), 
			   figs1h, labels = c(NA, "H"), ncol = 2, nrow = 1, widths = c(3,2.4), font.label = list(size = 20))

	test <- ggarrange(upper, lower, nrow = 2, ncol = 1, labels = c(NA, NA), heights = c(3,2.7))
	png(paste(sumDir, "figs1_cowplot.png", sep = "/"), res = figres, width = 36, height = 35, units = 'in')
	print(test)
	gar <- dev.off()

}

