## Get FlyBase Gene IDs used in am1
am1 <- read.table( "/home/intro-rna/2017/SRR891601.htseq", col.names = c( "id", "counts" ) )
head( am1 )


## Load biomaRt D. melanogaster dataset
library( "biomaRt" )
ensembl <- useMart( "ensembl" )
ensembl <- useDataset( "dmelanogaster_gene_ensembl", mart=ensembl )

# View( listAttributes( ensembl ) ) ... note trims after 1000 rows
listAttributes( ensembl )[ c(52,49,1,46,5,6:9,35:39,636,644), ]


## Retrieve look-up table
genenames <- getBM( attributes = c( "flybase_gene_id", "external_gene_name" ),
                    filter     = "flybase_gene_id", 
                    values     = am1$id, 
                    mart       = ensembl )
head( genenames )


## Annotate am1 data.frame
am1$name <- genenames$external_gene_name[ match( am1$id, genenames$flybase_gene_id ) ]
head( am1 )


## Find genes of interest
goi <- grepl( "PGRP", am1$name )
PGRPs <- am1[ goi, ]
PGRPs
