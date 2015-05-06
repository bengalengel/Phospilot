
suppressPackageStartupMessages(library("biomaRt"))
ensembl <- useMart(host = "grch37.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")

filenames <- list.files(path = "results/tables", pattern = "hg19.txt",
                        full.names = TRUE)

# Write header
cat("gene\ttranscript\tpeptide\tbiotype\teffect\tsnp\tref\talt\taf\tchr\tpos\tcdna\taa\n")

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
