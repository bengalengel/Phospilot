# Pipeline for cataloging Yoruba genetic variation

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

Versions used:

*  Python: 2.7
*  Java: 1.7.0_51
*  [SnpEff][]/[SnpSift][]: 4.1b
*  R: 3.1.1
*  [biomaRt][]: 2.22.0
*  Ensembl archive: 75 (uses GRCh37; http://feb2014.archive.ensembl.org)

[biomaRt]: http://bioconductor.org/packages/release/bioc/html/biomaRt.html
