#!/bin/bash

# Run IMPUTE2 result through SnpEff and keep coding sequence changes.

# Usage: bash run_snpeff.sh [ IMPUTE file ] [ Output directory ]
# Ex: bash run_snpeff.sh /mnt/lustre/data/internal/genotypes/hg19/YRI/chr21.hg19.impute2_haps.gz results/
# Note: Needs to be run from snpeff subdirectory so that it can find accessory scripts.

IMPUTE_FILE=$1
#IMPUTE_FILE=/mnt/lustre/data/internal/genotypes/hg19/YRI/chr21.hg19.impute2_haps.gz

OUTDIR=$2
mkdir -p $OUTDIR

# Software
PYTHON=/mnt/lustre/home/y.fourne/bin/virtualenv-1.11/py2.7/bin/python
CONVERT_VCF=./convert_impute2_to_vcf.py
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
# http://snpeff.sourceforge.net/SnpSift.html#filter
# See section "Effect prediction details" in http://snpeff.sourceforge.net/SnpEff_manual.html#input 
# for full list of potential effects.
# Use =~ operator to get regular expression matching in order to get joint effects like
# missense_variant&splice_region_variant and frameshift_variant&stop_gained.
mkdir -p $OUTDIR/snpeff_filter/
$JAVA -Xmx6g -jar $SNPSIFT filter \
  "(ANN[0].EFFECT =~ 'coding_sequence_variant') |
  (ANN[0].EFFECT =~ 'inframe_insertion') |
  (ANN[0].EFFECT =~ 'disruptive_inframe_insertion') |
  (ANN[0].EFFECT =~ 'inframe_deletion') |
  (ANN[0].EFFECT =~ 'disruptive_inframe_deletion') |
  (ANN[0].EFFECT =~ 'exon_loss_variant') |
  (ANN[0].EFFECT =~ 'frameshift_variant') |
  (ANN[0].EFFECT =~ 'gene_variant') |
  (ANN[0].EFFECT =~ 'missense_variant') |
  (ANN[0].EFFECT =~ 'initiatior_codon_variant') |
  (ANN[0].EFFECT =~ 'rare_amino_acid_variant') |
  (ANN[0].EFFECT =~ 'stop_lost') |
  (ANN[0].EFFECT =~ 'start_lost') |
  (ANN[0].EFFECT =~ 'stop_gained') |
  (ANN[0].EFFECT =~ 'transcript_variant')" \
  $OUTDIR/snpeff_oneperline/$BASE.snpeff.oneperline.vcf  \
  > $OUTDIR/snpeff_filter/$BASE.snpeff.filter.vcf

# Add functional predictions
# http://snpeff.sourceforge.net/SnpSift.html#dbNSFP
mkdir -p $OUTDIR/snpeff_dbnsfp/
$JAVA -Xmx6g -jar $SNPSIFT dbnsfp -v -a \
  -db dbNSFP.txt.gz \
  $OUTDIR/snpeff_filter/$BASE.snpeff.filter.vcf \
  > $OUTDIR/snpeff_dbnsfp/$BASE.snpeff.dbnsfp.vcf

# Extract fields
# http://snpeff.sourceforge.net/SnpSift.html#Extract
mkdir -p $OUTDIR/tables
$JAVA -Xmx6g -jar $SNPSIFT extractFields $OUTDIR/snpeff_dbnsfp/$BASE.snpeff.dbnsfp.vcf \
  "ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].EFFECT" "ID" "REF" "ALT" "AF" \
  "CHROM" "POS" "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
  "dbNSFP_Uniprot_acc" "dbNSFP_Interpro_domain" "dbNSFP_SIFT_pred" \
  "dbNSFP_Polyphen2_HDIV_pred" "dbNSFP_Polyphen2_HVAR_pred" \
  "dbNSFP_LRT_pred" "dbNSFP_MutationTaster_pred" "dbNSFP_phastCons100way_vertebrate" \
  "dbNSFP_1000Gp1_AF" "dbNSFP_1000Gp1_AFR_AF" "dbNSFP_1000Gp1_EUR_AF" "dbNSFP_1000Gp1_AMR_AF" \
  "dbNSFP_1000Gp1_ASN_AF" "dbNSFP_ESP6500_AA_AF" "dbNSFP_ESP6500_EA_AF"  "GEN[*].GT" \
  > $OUTDIR/tables/$BASE.txt
