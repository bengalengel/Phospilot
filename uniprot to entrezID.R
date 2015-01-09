##' Convert uniprot ID to entrez gene id
##'
##'
##' @title Uniprot2EG
##' @param uniprot uniprot id
##' @param organism only human supported now
##' @param getDescription get gene description or not
##' @param useBiomart using Biomart to query or using annotation packages.
##' @return a data.frame for mapping
##' @importFrom biomaRt getBM
##' @importFrom biomaRt useMart
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
Uniprot2EG <- function(uniprot, organism = "human", getDescription=FALSE, useBiomart=TRUE) {
  if (useBiomart == TRUE) {
    ##require(biomaRt)
    if (organism == "human") {
      wh_dataset = "hsapiens_gene_ensembl"
    } else if (organism == "mouse") {
      wh_dataset = "mmusculus_gene_ensembl"
    } else if (organism == "rat") {
      wh_dataset = "rnorvegicus_gene_ensembl"
    } else {
      stop ("Not supported yet...\n")
    }
    ensembl = useMart("ensembl", dataset=wh_dataset)
    if (getDescription) {
      ## UniProt/SwissProt Accession
      eg <- getBM(attributes=c("uniprot_swissprot_accession",
                               "entrezgene", "description"),
                  filters="uniprot_swissprot_accession",
                  values=uniprot, mart=ensembl)
      idx.notMap <- is.na(eg[,2])
      if( any(idx.notMap) ) {
        ## Uniprot/TrEMBL Accession
        ## uniprot_sptrembl
        eg2 <- getBM(attributes=c("uniprot_sptrembl",
                                  "entrezgene", "description"),
                     filters="uniprot_sptrembl",
                     values=as.character(eg[idx.notMap, 1]),
                     mart=ensembl)
      }
    } else {
      eg <- getBM(attributes=c("uniprot_swissprot_accession",
                               "entrezgene"),
                  filters="uniprot_swissprot_accession",
                  values=uniprot, mart=ensembl)
      idx.notMap <- is.na(eg[,2])
      if( any(idx.notMap) ) {
        ## Uniprot/TrEMBL Accession
        ## uniprot_sptrembl
        eg2 <- getBM(attributes=c("uniprot_sptrembl",
                                  "entrezgene", "description"),
                     filters="uniprot_sptrembl",
                     values=as.character(eg[idx.notMap, 1]),
                     mart=ensembl)
      }
    }
    names(eg)[1] <- "uniprot"
    if( any(idx.notMap) ) {
      names(eg2)[1] <- "uniprot"
      result <- rbind(eg, eg2)
    } else {
      result <- eg
    }
    result <- result[!is.na(result[,2]),]
    result <- result[result[,2] != "",]
    result <- unique(result)
  } else {
  }
  return (result)
}