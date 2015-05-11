# Confirm that snpeff does not change the reference and alternative allele
# from the original impute2 data.

# To submit from spudhead:
# for CHR in {1..22}
# do
#   echo "Rscript verify_alleles.R $CHR" | qsub -l h_vmem=8g -cwd -V -N verify.chr$CHR -o log -j y
# done

args <- commandArgs(trailingOnly = TRUE)
chr <- as.numeric(args[1])

print(chr)
# Read in snpeff table:
snpeff <- read.table(sprintf("results/tables/chr%d.hg19.txt", chr),
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# dim(snpeff)
colnames(snpeff)[1:12] <- c("gene", "transcript", "biotype", "effect", "snp",
                            "ref", "alt", "af", "chr", "pos", "cdna", "aa")
# snpeff[1, 1:13]
# Read in impute2 file:
impute2 <- read.table(gzfile(sprintf("/mnt/lustre/data/internal/genotypes/hg19/YRI/chr%d.hg19.impute2_haps.gz", chr)),
                      header = FALSE, sep = " ", stringsAsFactors = FALSE)
# dim(impute2)
colnames(impute2)[2:5] <- c("snp", "pos", "ref", "alt")
# impute2[1, 1:13]
# Verify all SNPs from snpeff in original impute2 file
stopifnot(length(intersect(snpeff$pos, impute2$pos)) ==
            length(unique(snpeff$pos)),
          length(intersect(snpeff$snp, impute2$snp)) ==
            length(unique(snpeff$snp[snpeff$snp != ""])))
# Clean:
# Remove the first column, which is just placeholder for chromosome
impute2_clean <- impute2[, -1]
# Convert unamed SNPs from "." to "" for consistency with snpeff
impute2_clean$snp <- ifelse(impute2$snp == ".", "", impute2$snp)
# Create unique id for merging
impute2_clean$id <- paste(impute2_clean$snp, impute2_clean$pos, sep = ".")
# Genotypes from snpeff are in format a1|a2. Need to separate into two
# adjacent columns.
snpeff_genos <- snpeff[, 13:ncol(snpeff)]
snpeff_allele_1 <- t(apply(snpeff_genos, 1, substr, start = 1, stop = 1))
snpeff_allele_2 <- t(apply(snpeff_genos, 1, substr, start = 3, stop = 3))
# Interleave the two alleles into one file where the first and second allele for
# an individual are adjacent columns.
snpeff_alleles <- matrix(nrow = nrow(snpeff_allele_1),
                         ncol = ncol(snpeff_allele_1) * 2)
for (i in 1:ncol(snpeff_allele_1)) {
  snpeff_alleles[, i * 2 - 1] <- snpeff_allele_1[, i]
  snpeff_alleles[, i * 2] <- snpeff_allele_2[, i]
}
snpeff_alleles <- apply(snpeff_alleles, 2, as.numeric)
snpeff_clean <- data.frame(snpeff[, c("snp", "pos", "ref", "alt")], snpeff_alleles,
                           stringsAsFactors = FALSE)
# Create unique id for merging
snpeff_clean$id <- paste(snpeff_clean$snp, snpeff_clean$pos, sep = ".")
# Merge the two files
merged <- merge(snpeff_clean, impute2_clean, by = "id",
                suffixes = c(".snpeff", ".impute2"))
stopifnot(nrow(merged) > 0)
# Confirm that the ref and alt alleles and the individual genotypes are
# identical.
for (i in 1:nrow(merged)) {
  if(!all(merged[i, 2:245] == merged[i, 246:489])) {
    if (interactive()) {
      browser()
    }
    cat(sprintf("Problem with SNP %s on chr%d\n", merged$snp.snpeff[i], chr))
  }
}
