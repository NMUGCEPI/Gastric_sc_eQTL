#!/bin/bash

type1=$1
celltype1=$2
type2=$3
celltype2=$4
Rscript /data/gc_sceqtl/sceqtl/eQTL/conditioningBetweenCellTypes_lm.R  ${type1} ${celltype1} ${type2} ${celltype2}
