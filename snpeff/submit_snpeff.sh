#!/bin/bash

mkdir -p log

for CHR in {1..22}
do
  echo "bash run_snpeff.sh /mnt/lustre/data/internal/Yoruba/IMPUTE/YRI/hg19/chr${CHR}.hg19.impute2_haps.gz results/" | qsub -l h_vmem=8g -cwd -V -N snpeff -o log -j y
done
