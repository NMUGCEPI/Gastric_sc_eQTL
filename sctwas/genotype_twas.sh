#!/bin/bash

type=$1
celltype=$2
plink --bfile /data/gc_sceqtl/twas/geno/ALLCHR_info05_SNP --noweb --keep /Public/wtp/blj/sceqtl/plinkfile208/${type}/keep_${celltype}id.sample --make-bed --out /data/gc_sceqtl/twas/geno/${type}/${celltype}_info05_SNP
plink --bfile /data/gc_sceqtl/twas/geno/${type}/${celltype}_info05_SNP --out /data/gc_sceqtl/twas/geno/${type}/${celltype}_info05_SNP_filter05  --maf 0.05 --hwe 0.001  --geno 0.05 --make-bed
