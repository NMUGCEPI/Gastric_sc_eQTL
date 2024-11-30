#!/bin/bash

type=$1
celltype=$2
plink --bfile /data/gc_sceqtl/sceqtl/plinkfile208/ALLCHR_info05 --noweb --keep /data/gc_sceqtl/sceqtl/plinkfile208/${type}/keep_${celltype}id.sample --make-bed --out /data/gc_sceqtl/herit/${celltype}/${celltype}_info05
plink --bfile /data/gc_sceqtl/herit/${celltype}/${celltype}_info05 --out /data/gc_sceqtl/herit/${celltype}/${celltype}_info05_filter05  --maf 0.05 --hwe 0.001  --geno 0.05 --make-bed
rm /data/gc_sceqtl/herit/${celltype}/${celltype}_info05.*
