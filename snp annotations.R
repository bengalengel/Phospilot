#here I am going to add SNP based annotation to multExpaned table for categorical enrichment test

##adding annotation to multExpanded1_withDE. The goal is to use a categorical enrichment test to see if phosphoproteins with nonsyn-snps unique to at least one of the lines are more likely to be over-represented in DE than those that do not have non-syn snps. 

# The next goal is to identify those sites with missing phosphorylated residues or regions surrounding the residues to get a sense of effect size after normalizing by protein.

#import yannick files
yannick18486 <- read.table("C:/Users/Brett/Desktop/pilot_18486_allNS.csv", sep = ',', header = T)
yannick18862 <- read.table("C:/Users/Brett/Desktop/pilot_18862_allNS.csv", sep = ',', header = T)
yannick19160 <- read.table("C:/Users/Brett/Desktop/pilot_19160_allNS.csv", sep = ',', header = T)
pqtl <- read.table("E:/My Documents/Dropbox/Postdoc-Gilad/Yannick SNP Phos/1260793_DatafileS1_pQTLs.csv", sep = ',', header = T)

#add a +/- based on UNIPROT.ID annotation, also add if isoform is present. Did I really add the annotation even if there is an
# isoform? NO see below
multExpanded1_withDE$Ind18486_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick18486$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick18862$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick19160$UNIPROT.ID, "+","-")

# test <- as.character(multExpanded1_withDE$Protein)
# isos <- multExpanded1_withDE[which(nchar(test)>6),]
# any(isos$Ind18486_SNP=="+")
# any(isos$Ind18486_SNPiso=="+")



#subsetting to first six characters to ensure isoforms get annoatated. Ignoring isoform maps from phosphopeptides for the moment
multExpanded1_withDE$Ind18486_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18486$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18862$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick19160$UNIPROT.ID, "+","-")


##subset to those phosphoobservations subject to DE (the background set of possibilities)
subtoDE <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",]
subtoDE <- subtoDE[!grepl("REV",subtoDE$Protein),]#now removed REVerse hits


##breakdown to the level of protein for categorical enrichment analysis
background <- subtoDE$Protein

#note that 5 of these proteins have a REVerse designation as the leading protein.
library(plyr)
count(grep("REV",subtoDE$Protein))
#they will be removed from this analysis for now, but I think I can go back and assign them to a different protein within the protein groups file later. This should also be done within the major file as well. I will ask about this in MQ help board.
background2 <- background[!grepl("REV", background)]
# > any(grep("REV",background2))
# [1] FALSE
# > any(grep("REV",background))
# [1] TRUE

#also note one protein is missing!
#missingprotein <- write.csv(background2,"background2.csv")


##**********************************using protien counts as background (not correcting for multiply phosphorylated proteins)
BGtable <- as.matrix(table(background2))#the parenthetical should be a factor vector of ids passed to this function
#remove entries with 0
BGtable <- BGtable[BGtable!=0,,drop=F]#1991 proteins

##add the snp annotation
BGtable <- as.data.frame(BGtable)
#now use the isoform information
anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPiso == "+" | multExpanded1_withDE$Ind18862_SNPiso == "+" | multExpanded1_withDE$Ind19160_SNPiso == "+",]
BGtable$anysnp <- ifelse(row.names(BGtable)%in%anysnp$Protein,"+","-")

# 
# anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNP == "+" | multExpanded1_withDE$Ind18862_SNP == "+" |                                  multExpanded1_withDE$Ind19160_SNP == "+",]
# BGtable$anysnp <- ifelse(row.names(BGtable)%in%anysnp$Protein,"+","-")

#somehow the first entry has no protein but is found with a SNP!

#also note the isoform issues. Correct annotation to map to any isoform which is a uniprot identifier followed by a dash.

# Now add the DE in any contrast information (already performed below)
BGtable$DEany <- ifelse(row.names(BGtable)%in%row.names(DEtable),"+","-")

#first row is DEany protein with and without a SNP
row1 <- c(nrow(BGtable[BGtable$DEany=="+" & BGtable$anysnp == "+",]), nrow(BGtable[BGtable$DEany == "+" & BGtable$anysnp == "-",])) 
#61 and 208 without isoforms
#96 and 173 with isoforms

#second row
row2 <- c(nrow(BGtable[BGtable$DEany=="-" & BGtable$anysnp == "+",]), nrow(BGtable[BGtable$DEany == "-" & BGtable$anysnp == "-",])) 
#290 and 1432 without isoforms
#481 and 1241 with isoforms

# This is going to be close...
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 0.006206
#two sided p = 0.01131


##****************************using full counts as background************************************8
anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPiso == "+" | multExpanded1_withDE$Ind18862_SNPiso == "+" | multExpanded1_withDE$Ind19160_SNPiso == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background2)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 
# 530 and 824 without specific isoform annotation (appropriate given this is bottom up MS) #note that there is a blank protein which may be causing the slight counting issue


#second row
row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 
# [1] 1072 2565

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 7.376e-11

##****************************using full counts as background and Yannicks full dataset***********8
#input Yannick's full dataset
#import yannick files
yannick18486full <- read.table("C:/Users/Brett/Desktop/pilot_18486_prot.csv",sep = ',', header = T)
yannick18862full <- read.table("C:/Users/Brett/Desktop/pilot_18862_prot.csv",sep = ',', header = T)
yannick19160full <- read.table("C:/Users/Brett/Desktop/pilot_19160_prot.csv",sep = ',', header = T)


#subsetting to first six characters to ensure isoforms get annoatated. Ignoring isoform maps from phosphopeptides for the moment
multExpanded1_withDE$Ind18486_SNPisofull = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18486full$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNPisofull = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18862full$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNPisofull = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick19160full$UNIPROT.ID, "+","-")

#add the phosphomotif information



anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPisofull == "+" | multExpanded1_withDE$Ind18862_SNPisofull == "+" | multExpanded1_withDE$Ind19160_SNPisofull == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background2)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 
# 530 and 824 without specific isoform annotation (appropriate given this is bottom up MS) #note that there is a blank protein which may be causing the slight counting issue
#now row1  558  796

#second row
row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 
# [1] 1072 2565
#now row2 1164 2473

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 9.664017e-10 

#some summary numbers at the site level
# contmatrix
             snp  nosnp
DE     row1  558  796
not DE row2 1164 2473
# Total sites that map to a protein with a non-syn snp and those that don't (first and second numbers respectively)
# > colSums(contmatrix)
# [1] 1722 3269
#34% of the phosphopeptides map to a protein that has a nonsyn snp

#contrast matrix for proteins
96 173
481 1241
protmatrix <- matrix(c(96,173,481,1241),nrow=2,ncol=2)
colSums(protmatrix)
[1]  269 1722
#13% of the proteins have a snp (yet 34% of the phosphopeptides map to a snp protein. proteins with snp tend to be longer and therefore have more phosphorylations)

##now for the tests of regions near the phosphosite***********************************************************888
#Are motif altering snps enriched for Diffphos?

anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPisofull == "+" | multExpanded1_withDE$Ind18862_SNPisofull == "+" | multExpanded1_withDE$Ind19160_SNPisofull == "+",]
#make background as data frame and add snp annotation
##breakdown to the level of protein for categorical enrichment analysis
background <- subtoDE$Protein
background <- as.data.frame(background)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#add the phosphovariant information to the anysnp dataframe, then subset this dataframe if variant present in any line
yannick18486motifvar <- yannick18486full[yannick18486full$Phosphomotif.positive..variant.is.within.7.Amino.acids.N.or.C.terminal.to.central.S.T.Y.=="yes",]##871 observations
yannick18862motifvar <- yannick18862full[yannick18862full$Phosphomotif.positive..variant.is.within.7.Amino.acids.N.or.C.terminal.to.central.S.T.Y.=="yes",]#487 observations
yannick19160motifvar <- yannick19160full[yannick19160full$Phosphomotif.positive..variant.is.within.7.Amino.acids.N.or.C.terminal.to.central.S.T.Y.=="yes",]#372 observations
###
anysnp$Ind18486_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18486motifvar$UNIPROT.ID, "+","-")
anysnp$Ind18862_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18862motifvar$UNIPROT.ID, "+","-")
anysnp$Ind19160_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick19160motifvar$UNIPROT.ID, "+","-")
#now subset to get all the 'motif'variants
motifvars <- anysnp[anysnp$Ind18486_motifproximal == "+" | anysnp$Ind18862_motifproximal == "+" | anysnp$Ind19160_motifproximal == "+",]
#add PROTEIN LEVEL proximal annotation to 'background'
background$nearmotif <- ifelse(background$Protein %in% motifvars$Protein, "+", "-")#isoforms will match but weren't important for mapping
## only 586 sites are annotated as near the motif out of 5000 mapping to 151 unique proteins
## add pqtl information
background$pqtl <- ifelse(substr(background$Protein,1,6) %in% pqtl$UNIPROT.ID, "+","-")## got 168/278!!

#enrichment of motif proximal variants at the protein level
#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$nearmotif == "+",]), nrow(background[background$DEany == "+" & background$nearmotif == "-",])) 
#second row
row2 <- c(nrow(background[background$DEany=="-" & background$nearmotif == "+",]), nrow(background[background$DEany == "-" & background$nearmotif == "-",])) 

# contmatrix
# [,1] [,2]
# row1  181 1173
# row2  405 3232

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
##just barely significant .01743108 

##Note that this assessment asks the question "are proteins with motif variants more likely than those without motif variants to be differentially phosphorylated?"

##"this is different than: "are phosphosites with a modification surrounding the residue more likely to be DiffPhos than those without a proximal modification?"

# The answers to these two questions speaks to the nature of PTM crosstalk.

#enrichment of pqtls in diff phos!!

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$pqtl == "+",]), nrow(background[background$DEany == "+" & background$pqtl == "-",])) 
#second row
row2 <- c(nrow(background[background$DEany=="-" & background$pqtl == "+",]), nrow(background[background$DEany == "-" & background$pqtl == "-",])) 

# [,1] [,2]
# row1   50 1304
# row2  118 3519

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
## Not significant! 0.2422782 this is just "is a protein that has a pqtl (that may not be different across individuals) enriched in DiffPhos"

##I need "are pqtls - where there is a difference in alleles between individuals in this subset - enriched in DE subset?"




#******************************************now for the protein normalized data using all snps and full counts**************************
##subset to those phosphoobservations subject to DE (the background set of possibilities)
subtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]#3603
subtoDEpn <- subtoDEpn[!grepl("REV",subtoDEpn$Protein),]#now removed REVerse hits and still 3603
##breakdown to the level of protein for categorical enrichment analysis
background <- subtoDEpn$Protein
#add the snp annotation for this analysis
anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPisofull == "+" | multExpanded1_withDE$Ind18862_SNPisofull == "+" | multExpanded1_withDE$Ind19160_SNPisofull == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDEpn$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDEpn[subtoDEpn$globalFsigpn=="+",]#972 obs with this truncated set
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#add the variant proximal annotation neglecting isoform information
anysnp$Ind18486_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18486motifvar$UNIPROT.ID, "+","-")
anysnp$Ind18862_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18862motifvar$UNIPROT.ID, "+","-")
anysnp$Ind19160_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick19160motifvar$UNIPROT.ID, "+","-")
#now subset to get all the 'motif'variants
motifvars <- anysnp[anysnp$Ind18486_motifproximal == "+" | anysnp$Ind18862_motifproximal == "+" | anysnp$Ind19160_motifproximal == "+",]

#add PROTEIN LEVEL proximity annotation to 'background'
background$nearmotif <- ifelse(background$Protein %in% motifvars$Protein, "+", "-")#isoforms will match but weren't important for mapping


#create the first row of the contingency matrix to test enrichment of snp containing proteins in DE subset.
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 

row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 


#FEtest
contmatrix <- rbind(row1,row2)
# [,1] [,2]
# row1  386  586
# row2  895 1736

result <- fisher.test(contmatrix, alternative = "g")
result$p.value
## 0.0009148777

##now for the tests of regions near the phosphosite
#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$nearmotif == "+",]), nrow(background[background$DEany == "+" & background$nearmotif == "-",])) 
#second row
row2 <- c(nrow(background[background$DEany=="-" & background$nearmotif == "+",]), nrow(background[background$DEany == "-" & background$nearmotif == "-",])) 
#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
##not significant this time.
## 0.4805825 

# > contmatrix
# [,1] [,2]
# row1  122  850
# row2  327 2304





