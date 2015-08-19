#vignette for customProDB
source("http://bioconductor.org/biocLite.R")
biocLite("customProDB")
library("customProDB")

browseVignettes("customProDB")

#an example of preparing the annotation using the UCSF genome browser and the package: "rtracklayer". Users must download manually the protein and coding sequence FASTA files. examles are loaded with the package
transcript_ids <- c("NM_001126112", "NM_033360", "NR_073499", "NM_004448", 
                    "NM_000179", "NR_029605", "NM_004333", "NM_001127511")
pepfasta <- system.file("extdata", "refseq_pro_seq.fasta", 
                        package="customProDB")
CDSfasta <- system.file("extdata", "refseq_coding_seq.fasta",
                        package="customProDB")
annotation_path <- tempdir()
PrepareAnnotationRefseq(genome='hg19', CDSfasta, pepfasta, annotation_path, 
                        dbsnp = NULL, transcript_ids=transcript_ids, 
                        splice_matrix=TRUE, COSMIC=FALSE)

#an example using the ensemble annotation and biomart's usemart

###################################################
### code chunk number 4: PrepareAnnoENSEMBL (eval = FALSE)
###################################################
## ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
##     host="feb2012.archive.ensembl.org", path="/biomart/martservice", 
##     archive=FALSE)
## annotation_path <- tempdir()
## transcript_ids <- c("ENST00000234420", "ENST00000269305", "ENST00000445888", 
##   		"ENST00000257430", "ENST00000457016", "ENST00000288602",
## 			"ENST00000269571", "ENST00000256078", "ENST00000384871")
## PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path, 
## 			splice_matrix=FALSE, dbsnp=NULL, 
##             transcript_ids=transcript_ids, COSMIC=FALSE)


# the next steps are used to build specific fasta files. 1) truncated main search using RPKM estimates from RNSseq to limit database
# 2) to use VCF files to produce a set of files specific to nonsyn and indels to add to the search 3) use the BED files to find existing and novel splice junctions. 

###################################################
### code chunk number 5: calculate
###################################################
load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))#exon annotation loaded
bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")#bam file
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))#
RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)

###################################################
### code chunk number 6: ouputpro
###################################################

#here they are outputting the FASTA
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
outf1 <- paste(tempdir(), '/test_rpkm.fasta', sep='')
Outputproseq(RPKM, 1, proteinseq, outf1, ids)


###################################################
### code chunk number 7: inpvcf
###################################################
# single sample
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)
length(vcf)
vcf[[1]][1:3]

###################################################
### code chunk number 8: inpvcfm
###################################################
# multiple samples in one VCF file
vcffile <- system.file("extdata", "test_mul.vcf", package="customProDB")
vcfs <- InputVcf(vcffile)

vcffile <- "D:/Yannick cluster data/SNPeffoutput files/Yoruba_filter.vcf"
testvcfs <- InputVcf(vcffile)


###################################################
### code chunk number 9: location
###################################################
table(values(vcf[[1]])[['INDEL']])
index <- which(values(vcf[[1]])[['INDEL']]==TRUE)
indelvcf <- vcf[[1]][index]

index <- which(values(vcf[[1]])[['INDEL']]==FALSE)
SNVvcf <- vcf[[1]][index]

load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite", package="customProDB"))
SNVloc <- Varlocation(SNVvcf,txdb,ids)
indelloc <- Varlocation(indelvcf,txdb,ids)
table(SNVloc[,'location'])


###################################################
### code chunk number 10: poscoding
###################################################
load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
load(system.file("extdata/refseq", "cosmic.RData", package="customProDB"))
postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding, COSMIC=cosmic)
postable_snv
postable_indel <- Positionincoding(indelvcf, exon)
postable_indel


###################################################
### code chunk number 11: aavar
###################################################
load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
txlist <- unique(postable_snv[, 'txid'])
codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist,]
mtab <- aaVariation (postable_snv, codingseq)
mtab


###################################################
### code chunk number 12: outvarpro
###################################################
outfile <-  paste(tempdir(), '/test_snv.fasta', sep='')
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
OutputVarproseq(mtab, proteinseq, outfile, ids)


###################################################
### code chunk number 13: outabber
###################################################
txlist_indel <- unique(postable_indel[, 'txid'])
codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] %in% txlist_indel, ]
outfile <-  paste(tempdir(), '/test_indel.fasta', sep='')
Outputaberrant(postable_indel, coding=codingseq_indel, proteinseq=proteinseq, 
               outfile=outfile, ids=ids)


###################################################
### code chunk number 14: junctype
###################################################
bedfile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")
jun <-  Bed2Range(bedfile,skip=1,covfilter=5)
jun

load(system.file("extdata/refseq", "splicemax.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
junction_type <- JunctionType(jun, splicemax, txdb, ids) 
junction_type[10:19,]
table(junction_type[, 'jun_type'])


###################################################
### code chunk number 15: novjunc
###################################################
outf_junc <- paste(tempdir(), '/test_junc.fasta',sep='')
library('BSgenome.Hsapiens.UCSC.hg19')
OutputNovelJun <- OutputNovelJun(junction_type, Hsapiens, outf_junc, 
                                 proteinseq)


###################################################
### code chunk number 16: sharedPrO
###################################################
path <- system.file("extdata/bams", package="customProDB")
bamFile<- paste(path, '/', list.files(path,pattern="*bam$"), sep='')
rpkms <- sapply(bamFile, function(x) 
  calculateRPKM(x, exon, proteincodingonly=TRUE, ids))
#colnames(rpkms) <- c('1', '2', '3')
#rpkms
outfile <- paste(tempdir(), '/test_rpkm_share.fasta', sep='')
pro <- OutputsharedPro(rpkms, cutoff=1, share_sample=2, proteinseq, 
                       outfile, ids)


###################################################
### code chunk number 17: mulvcf
###################################################
path <- system.file("extdata/vcfs", package="customProDB")
vcfFiles<- paste(path, '/', list.files(path, pattern="*vcf$"), sep='')
vcfs <- lapply(vcfFiles, function(x) InputVcf(x))
shared <- Multiple_VCF(vcfs, share_num=2)
shared


###################################################
### code chunk number 18: muljunc
###################################################
path <- system.file("extdata/beds", package="customProDB")
bedFiles<- paste(path, '/', list.files(path, pattern="*bed$"), sep='')
juncs <- lapply(bedFiles, function(x) Bed2Range(x, skip=1, covfilter=5))
sharedjun <- SharedJunc(juncs, share_num=2, ext_up=100, ext_down=100)
sharedjun


###################################################
### code chunk number 19: easyr
###################################################
bamFile <- system.file("extdata/bams", "test1_sort.bam", 
                       package="customProDB")
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
bedfile <- system.file("extdata", "junctions.bed", package="customProDB")
annotation_path <- system.file("extdata/refseq", package="customProDB")
outfile_path <- tempdir()
outfile_name='test'
easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, outfile_path, 
        outfile_name, rpkm_cutoff=1, INDEL=TRUE, lablersid=TRUE, COSMIC=TRUE, 
        nov_junction=FALSE) 


###################################################
### code chunk number 20: easyrmul
###################################################
bampath <- system.file("extdata/bams", package="customProDB")
vcfFile_path <- system.file("extdata/vcfs",  package="customProDB")
annotation_path <- system.file("extdata/refseq", package="customProDB")
outfile_path <- tempdir()	
outfile_name <- 'mult'
easyRun_mul(bampath, RPKM_mtx=NULL, vcfFile_path, annotation_path, rpkm_cutoff=1,
            share_num=2, var_shar_num=2, outfile_path, outfile_name, INDEL=TRUE,
            lablersid=TRUE, COSMIC=TRUE, nov_junction=FALSE)


###################################################
### code chunk number 21: FASTAnor
###################################################
outfile_path <- system.file("extdata/tmp", package="customProDB")
readLines(file(paste(outfile_path, '/test_rpkm.fasta', sep=''), 'rt'), 1)


###################################################
### code chunk number 22: FASTAsnv
###################################################
readLines(file(paste(outfile_path, '/test_snv.fasta', sep=''), 'rt'), 1)


###################################################
### code chunk number 23: FASTAindel
###################################################
readLines(file(paste(outfile_path, '/test_indel.fasta', sep=''), 'rt'), 1)


###################################################
### code chunk number 24: FASTAjun
###################################################
readLines(file(paste(outfile_path, '/test_junc.fasta', sep=''), 'rt'), 1)


###################################################
### code chunk number 25: SessionInfo
###################################################
sessionInfo()



