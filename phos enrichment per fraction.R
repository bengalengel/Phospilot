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
##summary is just like table
table(evidence$Modifications)

##make a global barplot (messy)
barplot(table(evidence$Modifications))


##make a binary barplot
x <- table(grepl("Phospho", evidence$Modifications))
rownames(x) <- c("non phosphorylated","phosphorylated")
x
barplot(x,ylab = "Total IDd peptides")

##now make a barplot cut by fraction; first a clue about fraction
class(evidence$Fraction)
allp <- table(evidence$Fraction)##I can make frequency table with integer as well as factor class variables
sum(evidence$Fraction)##doesn't add up
sum(allp)##adds up
##now the barplot
barplot(allp,ylab="all peptides IDd",xlab="fraction #")

##make a barplot of the phosphorylated peptides.

phos <- evidence[grepl("Phospho",evidence$Modifications),]
phospho <- table(phos$Fraction)
barplot(table(phos$Fraction), ylab = "phospho IDd", xlab = "fraction")


##make a barplot of each!
##first need to make a matrix out of two tables
parsed <- rbind(allp,phospho)

##side by sider
barplot(parsed,beside = T,ylab = "IDd peptides", xlab = "fraction", legend.text = T)
legend("topleft")

##stacked
nonphos <- allp-phospho##must subtract to get totals
parsedstk <- rbind(phospho,nonphos)##first table is plotted first
barplot(parsedstk,beside = F,ylab = "IDd peptides", xlab = "fraction", legend.text = T)

## there is a better way to do this at http://www.r-bloggers.com/stacked-bar-charts-in-r/

