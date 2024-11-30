##Cell composition GWAS
meta_all<-read.table('all_celltype_meta_info.xls',h=T)
dat <- table(meta_all$ID,meta_all$celltype)
dat <- as.matrix(as.data.frame.matrix(dat))
datp <- prop.table(dat,1)
datp <- as.data.frame(datp)
datp$ID <- rownames(datp)
#the arcsine square root transformation
datp$Arterial <- asin(sqrt(datp$Arterial))
datp$CD4 <- asin(sqrt(datp$CD4))
datp$CD8 <- asin(sqrt(datp$CD8))
datp$Chief <- asin(sqrt(datp$Chief))
datp$DCs <- asin(sqrt(datp$DCs))
datp$Endocrine <- asin(sqrt(datp$Endocrine))
datp$Enterocytes <- asin(sqrt(datp$Enterocytes))
datp$Inflam <- asin(sqrt(datp$Inflam))
datp$M1 <- asin(sqrt(datp$M1))
datp$M2 <- asin(sqrt(datp$M2))
datp$Mast <- asin(sqrt(datp$Mast))
datp$Myo <- asin(sqrt(datp$Myo))
datp$Naive <- asin(sqrt(datp$Naive))
datp$Neck <- asin(sqrt(datp$Neck))
datp$Parietal <- asin(sqrt(datp$Parietal))
datp$Pit <- asin(sqrt(datp$Pit))
datp$Plasma <- asin(sqrt(datp$Plasma))
datp$Stalk <- asin(sqrt(datp$Stalk))
datp$Tip <- asin(sqrt(datp$Tip))
write.table(datp,file='all_celltype_prop_log.xls',sep='\t',row.names=F,quote=F)

pheno <- fread('all_celltype_prop_log.xls')
fam <- fread('epi_info05_filter05.fam')
colnames(fam)[2] <- 'ID'
pheno <- merge(pheno,fam,by='ID')
pheno <- pheno[,c(21,1:20)]
colnames(pheno)[1:2] <- c('FID','IID') 
write.table(pheno,file='all_pheno.txt',sep='\t',row.names=F,col.names=T,quote=F)

###plink
plink --bfile epi_info05_filter05 \
      --pheno all_pheno.txt --all-pheno\
	  --covar covar.txt --allow-no-sex \
      --linear hide-covar --out celltype

#Skewed distribution converted into ordinal variables：Naive、Arterial、CD4、M1、M2、Enterocytes
Q20=quantile(pheno$Naive,seq(0.05,1,0.05))[[4]]
Q40=quantile(pheno$Naive,seq(0.05,1,0.05))[[8]]
Q60=quantile(pheno$Naive,seq(0.05,1,0.05))[[12]]
Q80=quantile(pheno$Naive,seq(0.05,1,0.05))[[16]]
pheno$Naive_20=0
pheno[pheno$Naive<=Q20,]$Naive_20=1
pheno[pheno$Naive>Q20&pheno$Naive<=Q40,]$Naive_20=2
pheno[pheno$Naive>Q40&pheno$Naive<=Q60,]$Naive_20=3
pheno[pheno$Naive>Q60&pheno$Naive<=Q80,]$Naive_20=4
pheno[pheno$Naive>Q80,]$Naive_20=5
Q20=quantile(pheno$Arterial,seq(0.05,1,0.05))[[4]]
Q40=quantile(pheno$Arterial,seq(0.05,1,0.05))[[8]]
Q60=quantile(pheno$Arterial,seq(0.05,1,0.05))[[12]]
Q80=quantile(pheno$Arterial,seq(0.05,1,0.05))[[16]]
pheno$Arterial_20=0
pheno[pheno$Arterial<=Q20,]$Arterial_20=1
pheno[pheno$Arterial>Q20&pheno$Arterial<=Q40,]$Arterial_20=2
pheno[pheno$Arterial>Q40&pheno$Arterial<=Q60,]$Arterial_20=3
pheno[pheno$Arterial>Q60&pheno$Arterial<=Q80,]$Arterial_20=4
pheno[pheno$Arterial>Q80,]$Arterial_20=5
Q20=quantile(pheno$CD4,seq(0.05,1,0.05))[[4]]
Q40=quantile(pheno$CD4,seq(0.05,1,0.05))[[8]]
Q60=quantile(pheno$CD4,seq(0.05,1,0.05))[[12]]
Q80=quantile(pheno$CD4,seq(0.05,1,0.05))[[16]]
pheno$CD4_20=0
pheno[pheno$CD4<=Q20,]$CD4_20=1
pheno[pheno$CD4>Q20&pheno$CD4<=Q40,]$CD4_20=2
pheno[pheno$CD4>Q40&pheno$CD4<=Q60,]$CD4_20=3
pheno[pheno$CD4>Q60&pheno$CD4<=Q80,]$CD4_20=4
pheno[pheno$CD4>Q80,]$CD4_20=5
Q20=quantile(pheno$M1,seq(0.05,1,0.05))[[4]]
Q40=quantile(pheno$M1,seq(0.05,1,0.05))[[8]]
Q60=quantile(pheno$M1,seq(0.05,1,0.05))[[12]]
Q80=quantile(pheno$M1,seq(0.05,1,0.05))[[16]]
pheno$M1_20=0
pheno[pheno$M1<=Q20,]$M1_20=1
pheno[pheno$M1>Q20&pheno$M1<=Q40,]$M1_20=2
pheno[pheno$M1>Q40&pheno$M1<=Q60,]$M1_20=3
pheno[pheno$M1>Q60&pheno$M1<=Q80,]$M1_20=4
pheno[pheno$M1>Q80,]$M1_20=5
Q20=quantile(pheno$M2,seq(0.05,1,0.05))[[4]]
Q40=quantile(pheno$M2,seq(0.05,1,0.05))[[8]]
Q60=quantile(pheno$M2,seq(0.05,1,0.05))[[12]]
Q80=quantile(pheno$M2,seq(0.05,1,0.05))[[16]]
pheno$M2_20=0
pheno[pheno$M2<=Q20,]$M2_20=1
pheno[pheno$M2>Q20&pheno$M2<=Q40,]$M2_20=2
pheno[pheno$M2>Q40&pheno$M2<=Q60,]$M2_20=3
pheno[pheno$M2>Q60&pheno$M2<=Q80,]$M2_20=4
pheno[pheno$M2>Q80,]$M2_20=5
Q20=quantile(pheno$Enterocytes,seq(0.05,1,0.05))[[4]]
Q40=quantile(pheno$Enterocytes,seq(0.05,1,0.05))[[8]]
Q60=quantile(pheno$Enterocytes,seq(0.05,1,0.05))[[12]]
Q80=quantile(pheno$Enterocytes,seq(0.05,1,0.05))[[16]]
pheno$Enterocytes_20=0
pheno[pheno$Enterocytes<=Q20,]$Enterocytes_20=1
pheno[pheno$Enterocytes>Q20&pheno$Enterocytes<=Q40,]$Enterocytes_20=2
pheno[pheno$Enterocytes>Q40&pheno$Enterocytes<=Q60,]$Enterocytes_20=3
pheno[pheno$Enterocytes>Q60&pheno$Enterocytes<=Q80,]$Enterocytes_20=4
pheno[pheno$Enterocytes>Q80,]$Enterocytes_20=5
pheno1 <- pheno[,c(1,2,22:33)]
write.table(pheno1,file='all_pheno_rank.txt',sep='\t',row.names=F,col.names=T,quote=F)

###plink
plink --bfile epi_info05_filter05 \
      --pheno all_pheno_rank.txt --all-pheno\
	  --covar covar.txt --allow-no-sex \
      --linear hide-covar --out celltype
