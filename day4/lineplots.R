selected_genes <- c(2:5)

library( "DESeq" )

# Use DESeq to normalize the read counts for each sample so they can be averaged
gut_metadata <- read.table( "/home/intro-rna/2017/midgut.tsv", header=TRUE )
cds <- newCountDataSetFromHTSeqCount( gut_metadata, directory="/home/intro-rna/2017" )
cds <- estimateSizeFactors( cds )

# make a data.frame of read counts in the correct orientation for the aggregate function
gene_counts <- data.frame(t(counts(cds, normalized=TRUE)[selected_genes,]))

# make a new data.frame of the mean count for each region 
mean_counts <- aggregate(gene_counts,by=list(gut_metadata$condition), mean)

# scale counts so they can be plotted on one graph
scaled_counts <- apply(mean_counts[,-1], 2, function(x) x/max(x))

plot(mean_counts[,2], type='l', xaxt='n', col="red", ylim=c(0,2000))
lines(mean_counts[,3], type='l', xaxt='n', col="green")
lines(mean_counts[,4], type='l', xaxt='n', col="blue")
lines(mean_counts[,5], type='l', xaxt='n', col="purple")
axis( 1, at=1:10, labels=mean_counts$Group.1)

plot(scaled_counts[,1], type='l', xaxt='n', col="red", ylim=c(0.0, 1.0))
lines(scaled_counts[,2], type='l', xaxt='n', col="green")
lines(scaled_counts[,3], type='l', xaxt='n', col="blue")
lines(scaled_counts[,4], type='l', xaxt='n', col="purple")
axis( 1, at=1:10, labels=mean_counts$Group.1)



