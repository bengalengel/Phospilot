#Script to calculate phospho IDs from single shot

#open the evidence file

## a bit abour filesystem management here. 
setwd("C:/Users//Brett/Documents/Pilot//3 x 2 x 2hr/txt/")


##the evidence file seems to be the one to use here as it contains all identified peptides 
evidence <- read.table("evidence.txt", sep = "\t", header=T, fill = T)

##it seems a little under half of my phospho data is enriched in total...
##how can a parse this into phospho vs unmodified/nonphospho
##summary is just like table
table(evidence$Modifications)

##make a global barplot (messy)
barplot(table(evidence$Modifications))


##make a binary barplot
x <- table(grepl("Phospho", evidence$Modifications))
rownames(x) <- c("non phosphorylated","phosphorylated")
x
barplot(x,ylab = "Total IDd peptides")

y <- as.matrix(x)
percent_enrichment <- y[2]/(y[1]+y[2])
