args = commandArgs(trailingOnly=TRUE)
cellLabel <- args[1]
print(cellLabel)
typeLabel <- args[2]
print(typeLabel)

wfilename <- sprintf('/data/gc_sceqtl/twas/weight/%s/',typeLabel)
bfilename <- sprintf('/data/gc_sceqtl/twas/geno/%s/%s_info05_SNP_filter05',cellLabel,typeLabel)
cfilename <- sprintf('/data/gc_sceqtl/twas/covar/%s_covar',typeLabel)

options(stringsAsFactors=F)
library(data.table)
library(R.utils)

region <- fread(sprintf("/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig/%s/%s/exp_%s_bed.bed.gz",cellLabel,typeLabel,typeLabel))
region <- region[,1:4]
colnames(region)[1] <- "Chr"

exp <- read.table(sprintf("/data/gc_sceqtl/twas/exp/%s_EXP_TWAS_gene_matrix",typeLabel))
exp <- cbind(ID = rownames(exp), exp)

region <- subset(region, ID %in% colnames(exp))
region$CHR <- region$Chr
region$flank_start <- region$start - 1e6
region$flank_start <- ifelse(region$flank_start<0, 0, region$flank_start)
region$flank_end <- region$end + 1e6

herit <- read.csv(sprintf("/data/gc_sceqtl/herit/%s/Heritability_sigresults_%s.csv",typeLabel,typeLabel))
herit <- subset(herit, Vg >=0 & Vp>=0 ) 
herit <- subset(herit, P < 0.1) 

## start calculation
for(i in 2:ncol(exp)){
	set.seed(1234)
	sub_range <- subset(region, ID==colnames(exp)[i])
	
	write.table(exp[,c(1,1,i)], paste0(wfilename, colnames(exp)[i], ".pheno"), col.names = F, row.names = F, quote = F, sep = "\t")
	commd <- paste0("plink --bfile ", bfilename, " --pheno ", wfilename, colnames(exp)[i], ".pheno --make-bed --out ", wfilename, colnames(exp)[i], " --chr ", sub_range$CHR, " --from-bp ", sub_range$flank_start, " --to-bp ", sub_range$flank_end)
	system(commd)
	rm(list = "commd")
	
	sub_herit <- subset(herit, gene==colnames(exp)[i])
	commd <- paste0("Rscript /data/blj/twas/fusion_twas-master/FUSION.compute_weights.R --bfile ", wfilename, colnames(exp)[i], " --tmp ", wfilename, colnames(exp)[i], ".tmp --covar ", cfilename, " --out ", wfilename, colnames(exp)[i], " --verbose 1 --save_hsq --PATH_gcta /data/blj/twas/fusion_twas-master/gcta_nr_robust --PATH_gemma gemma --PATH_plink plink --models top1,blup,lasso,enet --hsq_set ", sub_herit$Vg_p)
	system(commd)
	rm(list = "commd")
	
	commd <- paste0("rm ", wfilename, colnames(exp)[i], ".b*")
	system(commd)
	rm(list = "commd")

    commd <- paste0("rm ", wfilename, colnames(exp)[i], ".pheno")
	system(commd)
	rm(list = "commd")

    commd <- paste0("rm ", wfilename, colnames(exp)[i], ".log")
	system(commd)
	rm(list = "commd")

    commd <- paste0("rm ", wfilename, colnames(exp)[i], ".fam")
	system(commd)
	rm(list = "commd")
	
	print(i)
}

