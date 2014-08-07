##make a barchart of % of unique peptides identified in each fraction being phosphopeptides.

##the evidence file seems to be the one to use here as it contains all identified peptides 
evidence <- read.table("./MQ output/7_17_output/evidence.txt", sep = "\t", header=T, fill = T)
head(evidence)
names(evidence)
str(evidence$Fraction)
as.factor(evidence$Fraction)
str(evidence$Modifications)
summary(evidence$Modifications)
##it seems a little under half of my phospho data is enriched in total...
##how can a parse this into phospho vs unmodified/nonphospho

##now lets parse this by fraction

