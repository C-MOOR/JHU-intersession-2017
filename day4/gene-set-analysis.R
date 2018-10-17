## Identify differentially regulated genes
library( "DESeq" )
midgut <- read.table( "/home/intro-rna/2017/midgut.tsv", header=TRUE )
cds <- newCountDataSetFromHTSeqCount( midgut, directory="/home/intro-rna/2017" )
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
res <- nbinomTest( cds, "a1", "Fe" )


## Convert FlyBase Gene IDs to Entrez Gene IDs
library( "biomaRt" )
ensembl <- useMart( "ensembl" )
ensembl <- useDataset( "dmelanogaster_gene_ensembl", mart=ensembl )
conversion <- getBM( attributes = c( "flybase_gene_id", "external_gene_name", "entrezgene", "description" ),
                     filter     = "flybase_gene_id",
                     values     = res$id,
                     mart       = ensembl )

# Remove NA values and duplicates
a <- is.na( conversion$entrezgene )
b <- duplicated( conversion$flybase_gene_id ) | duplicated( conversion$flybase_gene_id, fromLast=TRUE )

res2 <- res[ res$id %in% conversion$flybase_gene_id[ !(a|b) ], ]
res2$entrez <- conversion$entrezgene[ match( res2$id, conversion$flybase_gene_id ) ]

# Remove six real/erroneous duplicate in biomaRt?
# res2[ duplicated( res2$entrez ), ]
# conversion[ grepl( "12798486", conversion$entrezgene ), ]
res2 <- res2[ !( duplicated( res2$entrez ) | duplicated( res2$entrez, fromLast=TRUE ) ), ]


## Gene set enrichment using goseq
library( "goseq" )

# Up-regulated genes
resUp <- res2$log2FoldChange > 0 & res2$padj < 0.10
names( resUp ) <- res2$entrez
resUp <- resUp[ !is.na( resUp ) ]

pwfUp <- nullp( resUp, "dm3", "refGene" )
GOup <- goseq( pwfUp, "dm3", "knownGene" )
GOup$padj <- p.adjust( GOup$over_represented_pvalue, method="BH" )

sigGOup <- GOup[ GOup$padj < 0.10, ]
summary( as.factor( sigGOup$ontology ) )
head( sigGOup[ sigGOup$ontology == "BP", 4:8 ] )

# Down-regulated genes
resDown <- res2$log2FoldChange < 0 & res2$padj < 0.10
names( resDown ) <- res2$entrez
resDown <- resDown[ !is.na( resDown ) ]

pwfDown <- nullp( resDown, "dm3", "refGene" )
GOdown <- goseq( pwfDown, "dm3", "knownGene" )
GOdown$padj <- p.adjust( GOdown$over_represented_pvalue, method="BH" )

sigGOdown <- GOdown[ GOdown$padj < 0.10, ]
summary( as.factor( sigGOdown$ontology ) )
head( sigGOdown[ sigGOdown$ontology == "BP", 6:8 ] )


