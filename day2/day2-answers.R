am4 <- read.table( "/home/intro-rna/SRR891604.htseq" )

dim( am4 )
head( am4 )
tail( am4 )

am4_a <- am4[ 1:17559, ]

totalreads <- sum( am4[,2] )
unaligned <- sum( am4[ 17560:17564, 2 ] )
unaligned / totalreads
summary( am4_a[,2] )

hist( log2( am4_a[,2] ), main="Reads mapping to genes in am4" )

sum( am4_a[,2] == 0 )

am30 <- read.table( "/home/intro-rna/SRR891630.htseq" )
am30_a <- am30[ 1:17559, ]
plot( am4_a[,2], am30_a[,2], log="xy" )
plot( am4_a[,2] + 1, am30_a[,2] + 1, log="xy", main="am4 vs am30 (log)" )



