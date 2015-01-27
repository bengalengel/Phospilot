#networKin relative enrichment barcharts
rm(list=ls(all=TRUE)) #start with empty workspace


# I think I should subset from the original networkin table in the future!!! Hard with the ID transfers...

#things to learn. readin multiple dataframe from folder....

#here background are the phosphosites subject to DE
NKbg <- read.table("NetKinout//networkin background.txt", sep = "\t", header=F, fill = T, quote = "")
#contrast specific reports
NKoutDE1 <- read.table("NetKinout//NKoutDE1.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE1down <- read.table("NetKinout//NKoutDE1down.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE1up <- read.table("NetKinout//NKoutDE1up.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE2 <- read.table("NetKinout//NKoutDE2.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE2down <- read.table("NetKinout//NKoutDE2down.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE2up <- read.table("NetKinout//NKoutDE2up.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE3 <- read.table("NetKinout//NKoutDE3.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE3down <- read.table("NetKinout//NKoutDE3down.txt", sep = "\t", header=F, fill = T, quote = "")
NKoutDE3up <- read.table("NetKinout//NKoutDE3up.txt", sep = "\t", header=F, fill = T, quote = "")



header <- c("substrate","position","id","networkin_score","tree","netphorest_group",
            "netphorest_score","string_identifier","string_score","substrate_name",
            "sequence","string_path")

colnames(NKbg) <- header
colnames(NKoutDE1) <- header
colnames(NKoutDE1down) <- header
colnames(NKoutDE1up) <- header
colnames(NKoutDE2) <- header
colnames(NKoutDE2up) <- header
colnames(NKoutDE2down) <- header
colnames(NKoutDE3) <- header
colnames(NKoutDE3up) <- header
colnames(NKoutDE3down) <- header


#Score truncation (1.5 TH) a likelihood score of 1 is neutral. Above 1 interactions are likely to be true.
NKBG <- NKbg[NKbg$networkin_score>=1.5,]
NKDE1 <- NKoutDE1[NKoutDE1$networkin_score>=1.5,]
NKDE1down <- NKoutDE1down[NKoutDE1down$networkin_score>=1.5,]
NKDE1up <- NKoutDE1up[NKoutDE1up$networkin_score>=1.5,]
NKDE2 <- NKoutDE2[NKoutDE2$networkin_score>=1.5,]
NKDE2down <- NKoutDE2down[NKoutDE2down$networkin_score>=1.5,]
NKDE2up <- NKoutDE2up[NKoutDE2up$networkin_score>=1.5,]
NKDE3 <- NKoutDE3[NKoutDE3$networkin_score>=1.5,]
NKDE3down <- NKoutDE3down[NKoutDE3down$networkin_score>=1.5,]
NKDE3up <- NKoutDE3up[NKoutDE3up$networkin_score>=1.5,]

#kinases
BGkinases <- NKBG[NKBG$tree=="KIN",]
DE1kinases <- NKDE1[NKDE1$tree=="KIN",]
DE1upkinases <- NKDE1up[NKDE1up$tree=="KIN",]
DE1downkinases <- NKDE1down[NKDE1down$tree=="KIN",]
DE2kinases <- NKDE2[NKDE2$tree=="KIN",]
DE2upkinases <- NKDE2up[NKDE2up$tree=="KIN",]
DE2downkinases <- NKDE2down[NKDE2down$tree=="KIN",]
DE3kinases <- NKDE3[NKDE3$tree=="KIN",]
DE3upkinases <- NKDE3up[NKDE3up$tree=="KIN",]
DE3downkinases <- NKDE3down[NKDE3down$tree=="KIN",]



#making contigency tables for each predicted kinase with examples



###########


NetEnrich <- function(x,y){
  #This function accepts lists of kinases found in a DE subset and background subset identified by the networkin program and return a list of adjusted pvalues for categorical enrichment using a fisher's exact test.x=background and y=enriched
  BGtable <- as.matrix(table(x))#the parenthetical should be a factor vector of ids passed to this function
  #remove entries with 0
  BGtable <- BGtable[BGtable!=0,,drop=F]#sum of table 2 is 4643
  DEtable <- as.matrix(table(y))#the parenthetical should be a passed DE factor vector of networKin output
  DEtable <- as.matrix(DEtable[row.names(DEtable) %in% row.names(BGtable),])#removing zeros and all factors not present in BG data
  #note the use of the %in% statement to control for rownames in DeTable but not in background
  #subset the background table in a similar way to ensure we are making the proper comparisons
  BGtable <- as.matrix(BGtable[row.names(BGtable) %in% row.names(DEtable),])
  NotDE <- BGtable-DEtable
  facttemp <- as.factor(row.names(DEtable))###########not quite working as a factor
  pvals <- c()
  #for each unique kinase in the DE1 kinases table I need to make a contingency table
  for(i in levels(facttemp)){#this would be the DE dataframe
    #make the first row of the contingency table
    if(DEtable[as.character(i),] & NotDE[as.character(i),] >= 0){
      DErow <- c(DEtable[as.character(i),],(sum(DEtable)-DEtable[as.character(i),]))
      NotDErow <- c(NotDE[as.character(i),],(sum(NotDE)-NotDE[as.character(i),]))
      contmatrix <- rbind(DErow,NotDErow)
      tmp <- fisher.test(contmatrix, alternative = "g")
      pvals <- c(pvals,tmp$p.value)
    }
  }
  ##multiple testing correction for pvals
  adjps <- p.adjust(pvals,method="BH")#none pass significance in this instance
  return(adjps)
}


##NetEnrich use
Background <- BGkinases$id
Enriched1 <- DE1kinases$id
Enriched2 <- DE2kinases$id
Enriched3 <- DE3kinases$id

Enriched1up <- DE1upkinases$id
Enriched1down <- DE1downkinases$id
Enriched2up <- DE2upkinases$id
Enriched2down <- DE2downkinases$id
Enriched3up <- DE3upkinases$id
Enriched3down <- DE3downkinases$id


#for testing
#x=Background
#y=Enriched2

DE1ps <- NetEnrich(Background,Enriched1)
DE2ps <- NetEnrich(Background,Enriched2)
DE3ps <- NetEnrich(Background,Enriched3)

DE1upps <- NetEnrich(Background,Enriched1up)
DE1downps <- NetEnrich(Background,Enriched1down)
DE2upps <- NetEnrich(Background,Enriched2up)
DE2downps <- NetEnrich(Background,Enriched2down)
DE3upps <- NetEnrich(Background,Enriched3up)
DE3downps <- NetEnrich(Background,Enriched3down)

#nothing enriched!!!! at alpha = .05 sig level for DE without direction contrasts

#Significance in the DE1down (5) contrast and the DE3down (1) contrast

> which(DE1downps<=.05)
[1] 19 36 42 46 70
> which(DE3downps<=.05)
[1] 60

#commonality is line #18486 I believe




#for each kinase I need a count and a fraction
barplot(table(as.character(DE1kinases$id)))
total <- sum(table(as.character(DE1kinases$id)))
DE1table <- as.matrix(table(as.character(DE1kinases$id)))
DE1table <- DE1table/total


myContingencyTable<-matrix(c(10,50,15,19925),nr=2)
fisher.test(myContingencyTable)

#I need to get the p-values then do multiple testing correction for each set with the entire proteome as background

summary(DE1kinases$id)
which(DE1kinases$id=="SOCS5")





