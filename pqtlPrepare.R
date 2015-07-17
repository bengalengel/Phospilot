
library(readxl)
library(R.utils)
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)

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
chain_file <- file.path(d, "hg18ToHg19.over.chain.gz")
if (!file.exists(chain_file)) {
  download.file(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz",
                destfile = chain_file)
}

pqtl_hg18 <- pqtl
colnames(pqtl_hg18)[6] <- "start"
pqtl_hg18$end <- pqtl_hg18$start
  
pqtl_hg18 <- makeGRangesFromDataFrame(pqtl_hg18,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE)

chain <- import.chain(gunzip(chain_file, overwrite = TRUE))

pqtl_hg19 <- liftOver(pqtl_hg18, chain)

# Obtain the genotypes from the VCF files
for (chr in 1:22) {
  vcf_name <- paste0("snpeff/results/vcf/chr", chr, ".hg19.vcf")
  file.exists(vcf_name)
}
