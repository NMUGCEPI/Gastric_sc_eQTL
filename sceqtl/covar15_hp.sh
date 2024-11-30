#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/sceqtl/eQTL/covar15_hp.r  ${type} ${celltype} 
