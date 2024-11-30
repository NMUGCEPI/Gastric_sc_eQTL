args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

#cellLabel <- "epi"
#typeLabel <- "pit"
options(stringsAsFactors=F)
library(data.table)
library(tidyfst)
library(janitor)
library(R.utils)

herit <- read.csv(sprintf("/data/gc_sceqtl/herit/%s/Heritability_sigresults_%s.csv",typeLabel,typeLabel)) 
sig_list <- herit %>%
    filter_dt((Vg >= 0) & (Vp >= 0) & (P < 0.1))
print(dim(sig_list))	
exp <- fread(sprintf("/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig/%s/%s/exp_%s_bed.bed.gz",cellLabel,typeLabel,typeLabel)) %>%
    filter_dt(ID %in% sig_list$gene) %>%
    t_dt() %>%
    row_to_names(row_number = 4)
write.table(exp, sprintf("/data/gc_sceqtl/twas/exp/%s_EXP_TWAS_gene_matrix",typeLabel), row.names = T, col.names = T, quote = F, sep = "\t") 
