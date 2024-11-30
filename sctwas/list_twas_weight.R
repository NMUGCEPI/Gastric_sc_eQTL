options(stringsAsFactors=F)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
number <- args[2]
print(number)

#cellLabel<-'chief'
#number<-'199'

list <- fread(sprintf("/data/gc_sceqtl/twas/allgene/list/%s_Exp_weights.list",cellLabel), header = F)
colnames(list)[1] <- "WGT"
list[, "ID" := tstrsplit(WGT, "/")[8]]
list[, "ID" := gsub(".wgt.RDat", "", ID)]

range <- fread("/data/gc_sceqtl/twas/Exp_region_for_heritability_2MB.list")
range <- range[match(list$ID, range$gene_id),]

weights <- merge(list, range, by.x = "ID", by.y = "gene_id")
weights$PANEL <- cellLabel
weights$WGT <- paste0(weights$ID, ".wgt.RDat")
weights$CHR <- weights$seqnames
weights$P0 <- weights$flank_start
weights$P1 <- weights$flank_end
weights$N <- number
weights <- weights[, c("PANEL", "WGT", "ID", "CHR", "P0", "P1", "N")]
write.table(weights, sprintf("/data/gc_sceqtl/twas/allgene/list/%s_Exp_weights_file.list",cellLabel), col.names = T, row.names = F, quote = F, sep = "\t")
