# Script for making a heatmap of the top 30 genes
# Make sure the plot window is big enough or it won't run

# Load data and prepare for clustering
library(DESeq)
gut_metadata <- read.table("/home/intro-rna/2017/midgut.tsv", header=TRUE)
cds <- newCountDataSetFromHTSeqCount(gut_metadata, directory="/home/intro-rna/2017")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
vsd <- varianceStabilizingTransformation(cds)

#Plot data
library(gplots)
# select top 30 genes
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
# plot heatmap of top 30 genes
heatmap.2(exprs(vsd)[select,], trace="none", margin=c(10, 6))
heatmap.2(counts(cds)[select,], trace="none", margin=c(10,6))

