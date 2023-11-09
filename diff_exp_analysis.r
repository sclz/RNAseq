# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0

# Librerie
library(DESeq2)
library(data.table)

library(ggplot2)
library(viridis)
library(hrbrthemes)

library(pheatmap)
library(tidyverse)

library_complexity_plot <- function(df) {
  tab = data.frame()  
  for (x in colnames(df)) {
    tab <- rbind(tab,c(x, "1", as.integer(sum(df[[x]] == 1))))
    tab <- rbind(tab,c(x, "2", as.integer(sum(df[[x]] == 2))))
    tab <- rbind(tab,c(x, "3", as.integer(sum(df[[x]] == 3))))
    tab <- rbind(tab,c(x, "5-10", as.integer(sum(df[[x]] >= 5 & df[[x]] <= 10))))
    tab <- rbind(tab,c(x, "10-20", as.integer(sum(df[[x]] > 10 & df[[x]] <= 20))))
    tab <- rbind(tab,c(x, "20-50", as.integer(sum(df[[x]] > 20 & df[[x]] <= 50))))
    tab <- rbind(tab,c(x, "50-100", as.integer(sum(df[[x]] > 50 & df[[x]] <= 100))))
    tab <- rbind(tab,c(x, ">100", as.integer(sum(df[[x]] > 100))))
  } 
  colnames(tab) <- c("Samples","Reads","Number of genes")
  tab$`Number of genes` <- as.integer(tab$`Number of genes`)
  ggplot(tab , aes(fill=tab$Reads, y=tab$`Number of genes`, x=tab$Samples)) + geom_bar(position="stack", stat="identity", orientation = "fill") + labs(x = "Samples") + labs(y = "Number of genes") + labs(fill = "Reads") + ggtitle("RNAseq QC") + scale_fill_viridis(discrete = T) + theme_ipsum()
  
}

# Ottieni gli argomenti dalla riga di comando
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript script.R count_table.tsv sample_info.tsv\n")
  q("no")
}

# Identificatore GEO
count_path <- args[1]
info_path <- args[2]
# Percorso del file dei conteggi
tbl <- as.matrix(fread(count_path, header = T, colClasses = "integer"), rownames = "GeneID")
#library_complexity_plot(as.data.frame(tbl))

# Percorso del file di annotazioni dei geni
apath <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
annot <- data.table::fread(apath, header = T, quote = "", stringsAsFactors = F, data.table = F)
rownames(annot) <- annot$GeneID

# Carica il file sample_info.tsv
sample_info <- read.table(info_path, header = T, stringsAsFactors = F)
#sample_info <- read.table("sample_info.tsv", header = T, stringsAsFactors = F)
sample_info$Group <- factor(sample_info$Group)

# filter out excluded samples
tbl <- tbl[, rownames(sample_info)]

# Pre-filtro dei geni a basso conteggio
keep <- rowSums(tbl >= 10) >= min(table(sample_info$Group))
tbl <- tbl[keep,]

# Crea il DESeqDataSet
ds <- DESeqDataSetFromMatrix(countData = tbl, colData = sample_info, design = ~Group)

ds <- estimateSizeFactors(ds) # Determine the size factors to use for normalization

normalized_counts <- counts(ds, normalized=TRUE) # Extract the normalized counts

#clustering e pca 

vsd <- vst(ds, blind=TRUE) # Transform the normalized counts

vsd_mat <- assay(vsd) # Extract the matrix of transformed counts

vsd_cor <- cor(vsd_mat) # Compute the correlation values between samples

pheatmap(vsd_cor, annotation = select(sample_info, Group)) # Plot the heatmap

plotPCA(vsd, intgroup="Group") # Plot the PCA of PC1 and PC2

# Esegui l'analisi di espressione genica differenziale
ds <- DESeq(ds, test = "Wald", sfType = "poscount")

# Estrai i risultati per i primi 250 geni
r <- results(ds, contrast = c("Group", levels(sample_info$Group)[1], levels(sample_info$Group)[2]), alpha = 0.05, pAdjustMethod = "fdr")
tT <- r[order(r$padj)[1:250],]

# Combina i risultati con le annotazioni dei geni
tT <- merge(as.data.frame(tT), annot, by = 0, sort = F)

# Seleziona colonne specifiche per l'output
tT <- subset(tT, select = c("GeneID", "padj", "pvalue", "lfcSE", "stat", "log2FoldChange", "baseMean", "Symbol", "Description"))

# Scrivi i risultati in un file tabulato
write.table(tT, file = "results.tsv", row.names = F, sep = "\t")

plotDispEsts(ds, main = "Dispersion Estimates")

# create histogram plot of p-values
hist(r$padj, breaks=seq(0, 1, length = 21), col = "grey", border = "white", 
     xlab = "", ylab = "", main = "Frequencies of padj-values")



# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2,1), cex.main=1.5)
plot(r$log2FoldChange, -log10(r$padj), main=paste(levels(sample_info$Group)[1], "vs", levels(sample_info$Group)[2]),
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)

# MD plot
par(mar=c(4,4,2,1), cex.main=1.5)
plot(log10(r$baseMean), r$log2FoldChange, main=paste(levels(sample_info$Group)[1], "vs", levels(sample_info$Group)[2]),
     xlab="log10(mean of normalized counts)", ylab="log2FoldChange", pch=20, cex=0.5)
with(subset(r, padj<0.05 & abs(log2FoldChange) >= 0),
     points(log10(baseMean), log2FoldChange, pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.05, sep=""), legend=c("down", "up"), pch=20,col=1:2)
abline(h=0)
palette(old.pal) # restore palette

################################################################
#   General expression data visualization
dat <- log10(counts(ds, normalized = T) + 1) # extract normalized counts

# box-and-whisker plot
lbl <- "log10(raw counts + 1)"
ord <- order(sample_info$Group)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
boxplot(dat[,ord], boxwex=0.6, notch=T, ylab="lg(norm.counts)", outline=F, las=2, col=sample_info$Group[ord])
legend("topleft", levels(sample_info$Group), fill=palette(), bty="n")

# UMAP plot (multi-dimensional scaling)
library(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
par(mar=c(3,3,2,6), xpd=TRUE, cex.main=1.5)
ump <- umap(t(dat), n_neighbors = 3, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=sample_info$Group, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(sample_info$Group), pch=20,
       col=1:length(levels(sample_info$Group)), title="Group", pt.cex=1.5)