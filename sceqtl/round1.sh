#!/bin/bash

type=$1
celltype=$2
Rscript /data/gc_sceqtl/sceqtl/eQTL/round1.R ${type} ${celltype}

