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
figs1_flag <- TRUE

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
	figs1j_file <- all_files[str_detect(all_files, "FigS1J")]

	figs1a <- readRDS(paste(sumDir, figs1a_file, sep = "/"))
	figs1b <- readRDS(paste(sumDir, figs1b_file, sep = "/"))
	figs1c <- readRDS(paste(sumDir, figs1c_file, sep = "/"))
	figs1d <- readRDS(paste(sumDir, figs1d_file, sep = "/"))
	figs1e <- readRDS(paste(sumDir, figs1e_file, sep = "/"))
	figs1f <- readRDS(paste(sumDir, figs1f_file, sep = "/"))
	figs1g <- readRDS(paste(sumDir, figs1g_file, sep = "/"))
	figs1h <- readRDS(paste(sumDir, figs1h_file, sep = "/"))
	figs1i <- readRDS(paste(sumDir, figs1i_file, sep = "/"))
	figs1j <- readRDS(paste(sumDir, figs1j_file, sep = "/"))


	figs1a <- figs1a + theme(text = element_text(size = 16))
	figs1b <- figs1b + theme(text = element_text(size = 16))
	figs1c <- figs1c + theme(text = element_text(size = 16))
	figs1d <- figs1d + theme(text = element_text(size = 16))
	figs1e <- figs1e + theme(text = element_text(size = 16))
	figs1f <- figs1f + theme(text = element_text(size = 16))
	figs1g <- figs1g + theme(text = element_text(size = 16))
	figs1h <- figs1h + theme(text = element_text(size = 16))
	figs1i <- figs1i + theme(text = element_text(size = 16))
	figs1j <- figs1j + theme(text = element_text(size = 16))

	upper <- ggarrange(
			   ggarrange(
				     ggarrange(figs1a, figs1b, figs1c, labels = c("A", "B", "C"), ncol = 3, nrow = 1, font.label = list(size = 20)), 
				     figs1e, labels = c(NA, "E"), ncol = 1, nrow = 2, heights = c(1,2.5), font.label = list(size = 20)), 
			   figs1d, labels = c(NA, "D"), ncol = 2, nrow = 1, widths = c(3,2.4), font.label = list(size = 20))

	lower <- ggarrange(
			   ggarrange(
				     ggarrange(figs1f, figs1g, figs1h, labels = c("F", "G", "H"), ncol = 3, nrow = 1, font.label = list(size = 20)), 
				     figs1j, labels = c(NA, "J"), ncol = 1, nrow = 2, heights = c(1,2.5), font.label = list(size = 20)), 
			   figs1i, labels = c(NA, "I"), ncol = 2, nrow = 1, widths = c(3,2.4), font.label = list(size = 20))

	test <- ggarrange(upper, lower, nrow = 2, ncol = 1, labels = c(NA, NA))
	png(paste(sumDir, "figs1_cowplot.png", sep = "/"), res = figres, width = 36, height = 36, units = 'in')
	print(test)
	gar <- dev.off()

}

