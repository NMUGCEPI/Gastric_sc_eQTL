#!/bin/bash
TYPE=$1
CHR=$2

echo ${TYPE}_${CHR}
Rscript /data/blj/twas/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /data/gc_sceqtl/twas/GC_gwas.sumstats \
    --weights /data/gc_sceqtl/twas/list/${TYPE}_Exp_weights_file.list \
    --weights_dir /data/gc_sceqtl/twas/weight/${TYPE} \
    --ref_ld_chr /data/blj/twas/LDREF_EAS/1000G.EAS. \
    --chr ${CHR} \
    --coloc_P 0.05 \
    --GWASN 21168 \
    --min_r2pred 0.4 \
    --out /data/gc_sceqtl/twas/result/${TYPE}/twas.chr${CHR}.dat | tee chr${CHR}.log