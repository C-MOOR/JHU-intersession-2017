metadata <- read.table( "/home/intro-rna/midgut.tsv", header=TRUE )

library( "DESeq" )

cds <- newCountDataSetFromHTSeqCount( metadata, directory="/home/intro-rna/" )
cds <- estimateSizeFactors( cds ) 
cds_t <- as.data.frame( t( counts( cds, normalized=TRUE ) ) ) 

library( "biomaRt" )

ensembl <- useMart( "ENSEMBL_MART_ENSEMBL" )
ensembl <- useDataset( "dmelanogaster_gene_ensembl", mart=ensembl )
genename <- getBM( attributes=c( "flybase_gene_id", "flybasename_gene" ),
                   filters="flybase_gene_id", values=colnames( cds_t ), mart=ensembl )
colnames( cds_t ) <- genename[ order( genename$flybase_gene_id ), ]$flybasename_gene

cds_t_a <- aggregate( cds_t, by=list( metadata$condition ), FUN=mean )
cds_final <- t( cds_t_a[ ,2:17560 ] )
colnames( cds_final ) <- cds_t_a[ ,1 ]
cds_final <- cds_final[ , unique( metadata$condition ) ]

plot(  cds_final[ "PGRP-SC1b", ], type='l', xaxt='n', col="green" )
lines( cds_final[ "PGRP-SC2" , ] * 0.025, xaxt='n', col="magenta" )
axis( at=1:10, side=1, labels=unique( metadata$condition ) )

write.csv( cds_final, "mean.counts.csv" )

