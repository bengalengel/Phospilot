
library(readxl)
library(R.utils)
library(GenomicRanges)
library(rtracklayer)

# Directory to save all data files produced by this script
d <- "."

# Download QTL results from Battle et al. 2014
xls <- file.path(d, "1260793_DatafileS1.xlsx")
if (!file.exists(xls)) {
  download.file(url = "http://www.sciencemag.org/content/suppl/2014/12/17/science.1260793.DC1/1260793_DatafileS1.xlsx",
                destfile = xls)
}

# Import pqtl data
pqtl <- read_excel(xls, sheet = 4, col_names = TRUE)

# Lift over from hg18 to hg19 coordinates
chain_file <- file.path(d, "hg18ToHg19.over.chain")
if (!file.exists(chain_file)) {
  chain_file_gz <- paste0(chain_file, ".gz")
  download.file(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz",
                destfile = chain_file_gz)
  gunzip(chain_file_gz)
}

pqtl_hg18 <- pqtl
colnames(pqtl_hg18)[6] <- "start"
pqtl_hg18$end <- pqtl_hg18$start

pqtl_hg18 <- makeGRangesFromDataFrame(pqtl_hg18,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE)

chain <- import.chain(chain_file)

pqtl_hg19 <- liftOver(pqtl_hg18, chain)
pqtl_hg19 <- as.data.frame(pqtl_hg19)

pqtl_final <- data.frame(ENSG = pqtl_hg19$ENSG,
                         chr = pqtl_hg19$seqnames,
                         perm.p.values = pqtl_hg19$perm.p.values,
                         snp.pvalue = pqtl_hg19$snp.pvalue,
                         snp.R.value = pqtl_hg19$snp.R.value,
                         hg19.pos = pqtl_hg19$start)
# Sort by chromosome (numeric) and position
pqtl_final <- pqtl_final[order(as.numeric(sub("chr", "", pqtl_final$chr)),
                               pqtl_final$hg19.pos), ]

head(pqtl_final)
write.table(pqtl_final, "pqtl-hg19.txt", quote = FALSE, sep = "\t", row.names = FALSE)
