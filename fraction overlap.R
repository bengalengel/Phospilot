
##count the number of peptides that are found in adjacent fractions.

peptides <- read.table("./MQ output/7_17_output/peptides.txt", sep = "\t", header=T, fill = T)


##some sort of nested for loop.
#for each observation
#for each slice
#if a value is present AND a value is present in the next slice
#add a value to count

##compare overrun into next fraction
##first must preallocate a matrix with expected iterations
iterations <- nrow(peptides)
variables <- 15

adjacency <- matrix(ncol = variables, nrow = iterations)
totunique <- matrix(ncol = variables, nrow = iterations)
j <- 1
for(i in 41:55){
  ##exists in fraction?
  totunique[,j] <- ifelse((!is.na(peptides[i])=="TRUE"),1,0)
    
#adjacency <- ifelse(!is.na(peptides[i]) && !is.na(peptides[i+1]),1,0)
  adjacency[,j] <- ifelse((!is.na(peptides[i]) & !is.na(peptides[i+1]) == "TRUE"),1,0)## winner winner. note && vs &!
  j <- j+1
}
##now turn adjacency and totunique into a data frame
adjacency <- data.frame(adjacency)
totunique <- data.frame(totunique)

##calculate sums on each fraction
olap <- unlist(lapply(adjacency,sum))
uniquefrac <- unlist(lapply(totunique,sum))

##plot total overlap by fraction
barplot(olap, xlab="fraction",ylab="overlapping peptides")

##plot of unique peptides per fraction
barplot(uniquefrac, xlab="fraction",ylab="unique peptides per fraction")##a bit unfair because a single experiment isn't being performd
##but may be indicative of what would happen with more populated samples

##stacked barplot of unique peptides found in successive fraction as a subset
uniquetofraction <- uniquefrac-olap##must subtract to get totals
parsedstk <- rbind(olap,uniquetofraction)##first table is plotted first
barplot(parsedstk,beside = F,ylab = "Unique Peptides Per Fraction", xlab = "Fraction", legend.text = T)

##barplot of percentage of each fraction found in the subsequent fraction
percentage <- olap/uniquefrac
barplot(percentage,ylab = "Percent Overlapping Per Fraction", xlab = "Fraction")

