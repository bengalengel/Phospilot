# Pipeline for cataloging Yoruba genetic variation

Download large data sets for annotation with [SnpEff][].

```bash
# Download genome
java -Xmx2g -jar snpeff download GRCh37.75
# Download dbNSFP for functional predictions
# dbNSFP.txt.gz is 9.1 GB!
# http://snpeff.sourceforge.net/SnpSift.html#dbNSFP
wget --no-check-certificate -O dbNSFP.txt.gz https://www.googledrive.com/host/0B7Ms5xMSFMYlSTY5dDJjcHVRZ3M
wget --no-check-certificate -O dbNSFP.txt.gz.tbi https://www.googledrive.com/host/0B7Ms5xMSFMYlOTV5RllpRjNHU2s
```

Start with YRI genotypes that have 1KG genotypes imputed with IMPUTE2.

Convert to VCF format, run through [SnpEff][] , and use [SnpSift][] to
only keep variants that affect coding sequence.

[SnpEff]: http://snpeff.sourceforge.net/index.html
[SnpSift]: http://snpeff.sourceforge.net/SnpSift.html

```bash
# Submit from head node
bash submit_snpeff.sh /mnt/lustre/data/internal/genotypes/hg19/YRI/ results/
```

This runs `run_snpeff.sh` for each chromosome file (only autosomes
available).

Combine the per-chromosome coding variants into one file, add Ensembl
peptide IDs, and add YRI samples names to header.

```bash
Rscript combine_snpeff.R results/tables /mnt/lustre/data/internal/genotypes/hg19/YRI/YRI_samples.txt > results/snpeff_final.txt
```

To verify that the SnpEff results designated the same alleles as the refernce and alternative versions as the original IMPUTE2 input data, run the following:

```bash
for CHR in {1..22}
do
 echo "Rscript verify_alleles.R $CHR" | qsub -l h_vmem=8g -cwd -V -N verify.chr$CHR -o log -j y
done
```

Check the log files for discordant SNPs:

```bash
cat log/verify*
```

Versions used:

*  Python: 2.7
*  Java: 1.7.0_51
*  [SnpEff][]/[SnpSift][]: 4.1b
*  R: 3.1.1
*  [biomaRt][]: 2.22.0
*  Ensembl archive: 75 (uses GRCh37; http://feb2014.archive.ensembl.org)

[biomaRt]: http://bioconductor.org/packages/release/bioc/html/biomaRt.html
