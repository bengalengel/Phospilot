ensembltest <- read.table("C:/Users//Brett/Desktop/mart_export.txt", sep = "\t", header = T)

any(duplicated(ensembltest))
duplicated(ensembltest$Ensembl.Transcript.ID)

dups <- which(duplicated(ensembltest$Ensembl.Transcript.ID))

duplicates <- ensembltest[dups,]

dups2 <- which(duplicated(ensembltest))

#all of the transcript duplicates are entry duplicates
intersect(dups,dups2)
setdiff(dups,dups2)

#all of the peptide duplicates are entry duplicates
dups3 <- which(duplicated(ensembltest$Ensembl.Protein.ID))
setdiff(dups2,dups3)


# All records..ehh
duplicates2 <- ensembltest[ensembltest$Ensembl.Transcript.ID == ensembltest$Ensembl.Transcript.ID[duplicated(ensembltest$Ensembl.Transcript.ID),],]
ALL_RECORDS <- df[df$ID==df$ID[duplicated(df$ID)],]
print(ALL_RECORDS) 