#!/usr/bin/env python

# Input the pQTL information. They are sorted by chromosome and position.
# This file was created by pqtlPrepare.R.
pqtl_file = open("pqtl-hg19.txt", "r")
pqtl_header = pqtl_file.readline()

# Input the names of the samples. We need the genotypes for 18486, 18862, and
# 19160.
samples = open("/mnt/lustre/data/internal/genotypes/hg19/YRI/YRI_samples.txt", "r")

samples.close()

out_file = open("pqtl-genos.txt", "w")
# Write header for output
out_file.write(pqtl_header.strip())

# Use a log file to investigate how the program is working
log_file = open("pqtl-log.txt", "w")

# Start searching on chr1
chr_current = "chr1"
vcf_file = open("snpeff/results/vcf/chr1.hg19.vcf", "r")
vcf_pos = None
log_file.write("current\t" + chr_current + "\n")

for pqtl in pqtl_file:
    pqtl_cols = pqtl.strip().split("\t")
    chr = pqtl_cols[1]
    pos = int(pqtl_cols[5])
    log_file.write(chr + "\t" + str(pos) + "\n")
    # If the pQTL is on the next chromosome, open its vcf file.
    if chr != chr_current:
        vcf_file.close()
        vcf_file = open("snpeff/results/vcf/" + chr + ".hg19.vcf", "r")
        chr_current = chr
        log_file.write("current\t" + chr_current + "\n")
    # If the pQTL is the same as the previous one, i.e. it is associated with
    # more than one gene, output the same SNP information instead of iterating
    # through more SNPs. This could also happen if the previous SNP could not be found and
    # the very next SNP that stopped the searching (because its position was
    # greater than that being searched for) was a pQTL.
    if pos == vcf_pos:
        out_file.write(pqtl.strip() + "\t" + vcf)
        continue
    # Search through the vcf file until the SNP is found. If the SNP position becomes
    # less than the current SNP in the vcf file, stop searching and output NA.
    for vcf in vcf_file:
        if vcf[0] == "#":
            continue
        vcf_cols = vcf.strip().split("\t")
        vcf_pos = int(vcf_cols[1])
        log_file.write(str(vcf_pos) + "\n")
        if pos == vcf_pos:
            out_file.write(pqtl.strip() + "\t" + vcf)
            break
        elif pos < vcf_pos:
            out_file.write(pqtl.strip() + "\t" + "NA" + "\n")
            break

pqtl_file.close()
vcf_file.close()
out_file.close()
log_file.close()
