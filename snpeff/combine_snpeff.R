# Combine per-chromosome snpeff results into one file.
# Add Ensembl peptide IDs.
# Add Yoruba sample IDs to header.

# Usage: Rscript combine_snpeff.R [ Input directory ] [ Samples file ]
# Ex: Rscript combine_snpeff.R results/tables /mnt/lustre/data/internal/genotypes/hg19/YRI/YRI_samples.txt > results/snpeff_final.txt
# Args:
#   Input directory: Tab-separated snpeff results from run_snpeff.sh
#   Samples file: Contains order of samples in IMPUTE2 files

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
samples_file <- args[2]

suppressPackageStartupMessages(library("biomaRt"))
ensembl <- useMart(host = "feb2014.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

filenames <- list.files(path = indir, pattern = "hg19.txt",
                        full.names = TRUE)

samples <- read.table(samples_file, header = FALSE, stringsAsFactors = FALSE)
# Test that all sample IDs start with NA
stopifnot(grepl("NA", samples[, 1]))

# Write header
cat("gene\ttranscript\tpeptide\tbiotype\teffect\tsnp\tref\talt\taf\tchr\tpos\tcdna\taa\tUniprot_acc\tInterpro_domain\tSIFT_pred\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_pred\tLRT_pred\tMutationTaster_pred\tphastCons100way_vertebrate\t1000Gp1_AF\t1000Gp1_AFR_AF\t1000Gp1_EUR_AF\t1000Gp1_AMR_AF\t1000Gp1_ASN_AF\tESP6500_AA_AF\tESP6500_EA_AF\t")
cat(samples[, 1], sep = "\t")
cat("\n")

for (f in filenames) {
  snps <- read.table(f, header = FALSE, sep = "\t", quote = "", skip = 1,
                     stringsAsFactors = FALSE)
  # Test that second column contains Ensembl transcript IDs
  stopifnot(grepl("ENST", snps[, 2]))
  
  # attributePages(ensembl)
  # atts <- listAttributes(ensembl, page = "feature_page")
  # atts[grep("Protein", atts$description), ]
  ensp <- getBM(attributes = c("ensembl_transcript_id", "ensembl_peptide_id"),
                filters = "ensembl_transcript_id",
                values = snps[, 2],
                mart = ensembl)
  # Merge the transcript and protein IDs
  snps_merged <- merge(snps, ensp, by.x = "V2", by.y = "ensembl_transcript_id",
                       all.x = TRUE)
  stopifnot(nrow(snps) == nrow(snps_merged))
  snps_final <- snps_merged[, c(2, 1, ncol(snps_merged), 3:ncol(snps))]
  stopifnot(grepl("ENSG", snps_final[, 1]),
            grepl("ENST", snps_final[, 2]),
            grepl("ENSP", snps_final[, 3]),
            ncol(snps_final) == ncol(snps) + 1)
  write.table(snps_final, file = "", quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}
