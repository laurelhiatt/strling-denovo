###need biocManager, biomaRt 
library(BiocManager)
library(biomaRt)
###library(EnsDb.Hsapiens.v79) This is the option below; recommended on stackoverflow as most accurate
### If interested, separate syntax is necessary.
getwd()
setwd("/Users/quinlan/Documents/Git/strling-denovo") ###or wherever you want to be
list.files() ###make sure the file you want is where you want it
data <- read.csv("STRs150bpknowngenesoverlapSORTEDSIGUNIQUE", stringsAsFactors = TRUE, header = FALSE, sep='\t')
ensembl.genes <- data$V19 ##filter by whatever row has your ensembl info
print(ensembl.genes) ###make sure you got what you wanted

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

G_list <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_gene_id","hgnc_symbol"),values= ensembl.genes, mart=mart)
### The filter must correspond with ensemble.genes, and the attributes determine our G_list output
### We can merge here if we would like, dependent on our data dataframe and goals.

write.table(G_list, file='nameoffile.tsv', quote=FALSE, index = FALSE, sep='\t')

