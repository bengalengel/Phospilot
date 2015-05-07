#!/bin/bash

# Submit one job per chromosome to run with run_snpeff.sh.

# Usage: bash submit_snpeff.sh [ Input directory ] [ Output directory ]
# Ex: bash submit_snpeff.sh /mnt/lustre/data/internal/genotypes/hg19/YRI/ results/
# Note: Needs to be run from snpeff subdirectory so that it can find accessory scripts.
# Note: log files will be written to the subdirectory log/ in the current working directory
# Args:
#  Input directory: Contains IMPUTE2 output files with filenames like chr21.hg19.impute2_haps.gz
#  Output directory: Directory to write results.

INDIR=$1
OUTDIR=$2

mkdir -p log

for CHR in {1..22}
do
  echo "bash run_snpeff.sh $INDIR/chr${CHR}.hg19.impute2_haps.gz $OUTDIR" | qsub -l h_vmem=8g -cwd -V -N chr.${CHR}.snpeff -o log -j y
done
