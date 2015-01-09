#networKin relative enrichment barcharts



#things to learn. readin multiple dataframe from folder....
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

colnames(NKoutDE1) <- header
colnames(NKoutDE1down) <- header
colnames(NKoutDE1up) <- header
colnames(NKoutDE2) <- header
colnames(NKoutDE2up) <- header
colnames(NKoutDE2down) <- header
colnames(NKoutDE3) <- header
colnames(NKoutDE3up) <- header
colnames(NKoutDE3down) <- header


#score truncation (1.5 TH)
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
DE1kinases <- NKDE1[NKDE1$tree=="KIN",]
DE1upkinases <- NKDE1up[NKDE1up$tree=="KIN",]
DE1downkinases <- NKDE1down[NKDE1down$tree=="KIN",]
DE2kinases <- NKDE2[NKDE2$tree=="KIN",]
DE2upkinases <- NKDE2up[NKDE2up$tree=="KIN",]
DE2downkinases <- NKDE2down[NKDE2down$tree=="KIN",]
DE3kinases <- NKDE3[NKDE3$tree=="KIN",]
DE3upkinases <- NKDE3up[NKDE3up$tree=="KIN",]
DE3downkinases <- NKDE3down[NKDE3down$tree=="KIN",]

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





