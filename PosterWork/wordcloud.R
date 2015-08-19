install.packages(c("wordcloud","tm"),repos="http://cran.r-project.org")
library(wordcloud)
library(tm)



wordcloud("This is it now dude")

names <- c("hat","mat","couch","douche","friendly")

wordcloud(names, freq = c(8,2,4,6,8.8), min.freq = 1, random.order=F, colors="red")

perDE1 <- read.table("perseus enrichment/DEcont1.txt", sep = "\t", header=T, fill = T, quote = "")
perDE1down <- read.table("perseus enrichment/cont1down.txt", sep = "\t", header=T, fill = T, quote = "")
perDE1up <- read.table("perseus enrichment/cont1upenrich.txt", sep = "\t", header=T, fill = T, quote = "")
perDE2 <- read.table("perseus enrichment/DEcont2.txt", sep = "\t", header=T, fill = T, quote = "")
perDE2down <- read.table("perseus enrichment/cont2down.txt", sep = "\t", header=T, fill = T, quote = "")
perDE2up <- read.table("perseus enrichment/cont2up.txt", sep = "\t", header=T, fill = T, quote = "")
perDE3 <- read.table("perseus enrichment/DEcont3.txt", sep = "\t", header=T, fill = T, quote = "")
perDE3down <- read.table("perseus enrichment/cont3down.txt", sep = "\t", header=T, fill = T, quote = "")
perDE3up <- read.table("perseus enrichment/cont3up.txt", sep = "\t", header=T, fill = T, quote = "")

#Subset
motifsDE1 <- perDE1[perDE1$Category.column=="Motifs" & perDE1$Selection.value == "+",]
motifsDE2 <- perDE2[perDE2$Category.column=="Motifs" & perDE2$Selection.value == "+",]
motifsDE3 <- perDE3[perDE3$Category.column=="Motifs" & perDE3$Selection.value == "+",]

names <- as.character(motifsDE3$Category.value)
names
#remove "substrate binding & motif"
library(stringr)
names <- gsub("binding","",names)
names <- gsub("substrate","",names)
names <- gsub("motif","",names)
names  <- str_trim(names)#remove whitespace at the end

p <- -log10(motifsDE3$P.value)


wordcloud(names, freq = p2, min.freq = 1, random.order=F, colors=brewer.pal(8, "Dark2"))



