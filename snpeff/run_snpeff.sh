#!/bin/bash

# Usage: bash run_snpeff.sh [ IMPUTE file ] [ Output directory ]
# Ex: bash run_snpeff.sh /mnt/lustre/data/internal/Yoruba/IMPUTE/YRI/hg19/chr21.hg19.impute2_haps.gz results/

IMPUTE_FILE=$1
#IMPUTE_FILE=/mnt/lustre/data/internal/Yoruba/IMPUTE/YRI/hg19/chr21.hg19.impute2_haps.gz

OUTDIR=$2
mkdir -p $OUTDIR

# Software
PYTHON=/mnt/lustre/home/y.fourne/bin/virtualenv-1.11/py2.7/bin/python
CONVERT_VCF=/mnt/gluster/home/jdblischak/phospho/convert_impute2_to_vcf.py
JAVA=/mnt/lustre/home/jdblischak/src/java/jre1.7.0_51/bin/java
SNPEFF=/mnt/lustre/home/jdblischak/src/snpEff/snpEff.jar
ONE_EFFECT_PER_LINE=/mnt/lustre/home/jdblischak/src/snpEff/scripts/vcfEffOnePerLine.pl
SNPSIFT=/mnt/lustre/home/jdblischak/src/snpEff/SnpSift.jar

BASE=`basename $IMPUTE_FILE`
BASE=${BASE%.impute2_haps.gz}
# echo $BASE

# Convert from impute2 to vcf format
mkdir -p $OUTDIR/vcf
$PYTHON $CONVERT_VCF $IMPUTE_FILE > $OUTDIR/vcf/$BASE.vcf

# Run snpEff
mkdir -p $OUTDIR/snpeff
$JAVA -Xmx6g -jar $SNPEFF -v \
  -no-downstream -no-intergenic -no-intron -no-upstream -no-utr \
  GRCh37.75 $OUTDIR/vcf/$BASE.vcf > $OUTDIR/snpeff/$BASE.snpeff.vcf

# Good idea to download databases beforehand
#JAVA -Xmx2g -jar $SNPEFF download GRCh37.75
#JAVA -Xmx2g -jar $SNPEFF download hg19

# Reduce to one effect per line
mkdir -p $OUTDIR/snpeff_oneperline
cat $OUTDIR/snpeff/$BASE.snpeff.vcf | $ONE_EFFECT_PER_LINE  > $OUTDIR/snpeff_oneperline/$BASE.snpeff.oneperline.vcf

# Filter
mkdir -p $OUTDIR/snpeff_filter/
$JAVA -Xmx6g -jar $SNPSIFT filter \
"(ANN[0].EFFECT == 'missense_variant') | (ANN[0].EFFECT == 'disruptive_inframe_deletion') | (ANN[0].EFFECT == 'frameshift_variant') | (ANN[0].EFFECT == 'frameshift_variant&stop_gained') | (ANN[0].EFFECT == 'inframe_deletion') | (ANN[0].EFFECT == 'inframe_insertion') | (ANN[0].EFFECT == 'missense_variant&splice_region_variant') | (ANN[0].EFFECT == 'start_gained') | (ANN[0].EFFECT == 'start_lost') | (ANN[0].EFFECT == 'stop_gained') | (ANN[0].EFFECT == 'stop_lost')" \
$OUTDIR/snpeff_oneperline/$BASE.snpeff.oneperline.vcf  \
> $OUTDIR/snpeff_filter/$BASE.snpeff.filter.vcf

# Extract fields
mkdir -p $OUTDIR/tables
$JAVA -Xmx6g -jar $SNPSIFT extractFields $OUTDIR/snpeff_filter/$BASE.snpeff.filter.vcf \
"ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ID" "REF" "ALT" "AF" "CHROM" "POS" "ANN[*].HGVS_P" "GEN[*].GT" \
> $OUTDIR/tables/$BASE.txt
