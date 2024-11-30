# Call the cell type
args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

library(data.table)

# Directory paths
output.dir <- "/data/gc_sceqtl/sceqtl/eQTL"

# Input filenames
input <- sprintf("%s/%s/%s_hp/covar/%s.combined_covariates.txt",output.dir,cellLabel,typeLabel,typeLabel)
output <- sprintf("%s/%s/%s_hp/covar/%s.combined_covariates15.txt",output.dir,cellLabel,typeLabel,typeLabel)

covar<-fread(input)
write.table(covar[c(1:12,26,27,28),],output,col.names=T,row.names=F,quote=F,sep='\t')

q()
