#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/sceqtl/eQTL/celltype_sig_eqtl.R  ${type} ${celltype} 
