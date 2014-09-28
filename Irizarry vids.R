# Rafa Irizarry notes from statistics for genomics vids

# log transform the data
data <- log2(multExpanded[,newnames])

# Rafael has replaced the row names with a useful identifier to play with this dataset on its own. Perhaps this is a worthwhile thing to do for me
# as well

#look at the histograms of the data note some pretty low values here
hist(data[,2], breaks=16)

#summarize the data
summary(data[,2])

# can summarize the entire data frame
summary(data)

# Can do the row means to get average per biological samples
HL16770 <- rowMeans(data[,1:2])
HL16778 <- rowMeans(data[,3:4])
HL16788 <- rowMeans(data[,5:6])

##MA plot derived from first two samples (ratio on y and average value of protein on A)
A <- rowMeans(data[,1:4])
plot(A,HL16778-HL16770)

##which protein has the biggest differences
which.min(HL16778-HL16770)##gives row name and position within vector

# get the top 10

o <- order(HL16778-HL16770) ##can use abs for absolute vaue and  change flag for decreasing to get biggest changes in either direction.
rownames(data)[o[1:20]]##the order vector returns the position within the original vector and rownames returns the name for that row within the 'data' dataframe

# Is the difference big due to actually change or does the gene simply vary alot across samples?
# Calculate variance

# First must install genefilter
source("http://bioconductor.org/biocLite.R")##I have run a script here using 'source' Note how to get new version when I update R?BiocUpgrade
biocLite("genefilter")

# installed
# package ‘BiocGenerics’ successfully unpacked and MD5 sums checked
# package ‘DBI’ successfully unpacked and MD5 sums checked
# package ‘RSQLite’ successfully unpacked and MD5 sums checked
# package ‘IRanges’ successfully unpacked and MD5 sums checked
# package ‘xtable’ successfully unpacked and MD5 sums checked
# package ‘XML’ successfully unpacked and MD5 sums checked
# package ‘AnnotationDbi’ successfully unpacked and MD5 sums checked
# package ‘annotate’ successfully unpacked and MD5 sums checked
# package ‘Biobase’ successfully unpacked and MD5 sums checked
# package ‘genefilter’ successfully unpacked and MD5 sums checked


library(genefilter)
# Get the sd for each gene
s1 <- rowSds(data[,1:2])
s2 <- rowSds(data[,3:4])

# Calculating a t-test manually... note that there is a function rowttest! Note my sample n is two replicates/condition.
ttest <- (HL16778-HL16770)/sqrt(s1^2/2 + s2^2/2)

ttest2 <- (HL16770-HL16778)/sqrt(s1^2/2 + s2^2/2)

# ttest2 <- rowttests(data)
hist(ttest, breaks = 100)#note there are some crazy large differences here
hist(ttest2, breaks = 100)

n <- ttest[!is.na(ttest)]#note this is for a vector

#Calculating pvalues ####################################################
pval <- 2*(1-pt(abs(ttest),2))   #pt Returns the chance that X>x; for two sided must multiply by 2!; df is # samples in ea condition -2; 

pt returns the proportion of the tdist that is greater than the imputted value; for example
pt(0,2) #returns 1/2 
pt(5,2) #returns a large number

# Here I have used the abs value so that the negative tstats return the correct value. I am multiplying by 2 to make it a two sided

pval2 <- 2*(1-pt(abs(ttest2),2))  
pval3 <- 2*(1-pt(abs(n),2))  #check for na issues
##########################################################################################


###########################QC of p-values######################## 

hist(pval)##shit there are many DE phospho?!!
hist(pval2)##these are identical of course (just checking...)
hist(pval3)#note that hist removes the NAs

# I should see a relatively flat distribution here!!!


# QQplots of t statistic for QC
qqt(n,df=2)
abline(0,1)











# Volcano plot
plot(HL16778-HL16770, -log10(pval), main="volcano")





