## Load biomaRt D. melanogaster dataset
library( "biomaRt" )
ensembl <- useMart( "ensembl" )
ensembl <- useDataset( "dmelanogaster_gene_ensembl", mart=ensembl )


## Get FlyBase Gene IDs used in am1
am1 <- read.table( "/home/intro-rna/2017/SRR891601.htseq" )


## Retrieve all Gene Ontology terms used to tag a gene in am1
goterms <- getBM( attributes = c( "go_id", "name_1006" ),
                  filter     = "flybase_gene_id",
                  values     = am1[ , 1 ],
                  mart       = ensembl )
head( goterms )


## Find go_id for a gene set of interest
goi <- grepl( "primary miRNA", goterms$name_1006 )
goterms[ goi, ]


## Find genes tagged with go_id 
genes <- getBM( attributes = c( "external_gene_name", "name_1006", "description" ),
                filter     = "go_id",
                values     = goterms[ goi, 1 ],
                mart       = ensembl )
head( genes )

