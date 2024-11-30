# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

# Example
# cellLabel <- "epi"
# typeLabel <- "chief"

library(data.table)
library(qvalue)
library(tidyverse)
library(broom)
library(dplyr)

# Directory paths
output.dir <- "/data/gc_sceqtl/sceqtl/eQTL"

# Input filenames
qtl_filename <- sprintf("%s/%s/%s_hp/%s_eqtl.allpairs.txt", output.dir, cellLabel, typeLabel, typeLabel)

gtf_gene<-read.csv('/data/gc_sceqtl/sceqtl/eQTL/GTF/gtf_gene_auto.csv',header=T)

qtl = fread(qtl_filename)
pvalues <- qtl$pval_nominal
qobj <- qvalue(p = pvalues, pi0=1)
qtl <- qtl %>% add_column(qvalue=qobj$qvalues,localFDR=qobj$lfdr)
qtl_sig <- qtl[(qtl$qvalue < 0.05),]
qtl_sig<-merge(qtl_sig,gtf_gene,by='gene_id')
qtl_sig1<- qtl_sig[(qtl_sig$qvalue < 0.01),]
print(length(unique(qtl_sig$gene_name)))
print(length(unique(qtl_sig1$gene_name)))
write.table(qtl_sig,sprintf("%s/%s/%s_hp/%s_qtl_sig_FDR05.xls", output.dir, cellLabel, typeLabel, typeLabel),quote = F,sep = "\t",row.names = F,col.names = T)

print('JOB IS DONE!')

quit()


