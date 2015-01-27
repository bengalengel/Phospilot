#here I am going to add SNP based annotation to multExpaned table for categorical enrichment test

##adding annotation to multExpanded1_withDE. The goal is to use a categorical enrichment test to see if proteins with snps are more likely to be over-represented in DE than those that do not have non-syn snps. 

# The next goal is to identify those sites with missing phosphorylated residues or regions surrounding the residues to get a sense of effect size after normalizing by protein.

#add a +/- based on 18486 UNIPROT.ID annotation, also add if isoform is present
multExpanded1_withDE$Ind18486_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick18486$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick18862$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick19160$UNIPROT.ID, "+","-")


#limit to the first six characters because after this is a dash and an isoform (ignoring REV enrtries for the moment)
multExpanded1_withDE$Ind18486_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18486$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18862$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick19160$UNIPROT.ID, "+","-")


##subset to those phosphoobservations subject to DE (the background set of possibilities)
subtoDE <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",]

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





# Now for the DE in any contrast vector of protein names
DEany <- multExpanded1_withDE[multExpanded1_withDE$DEcont1=="+"|multExpanded1_withDE$DEcont2=="+"|multExpanded1_withDE$DEcont3=="+",]
#remove the reverse for now
DEany <- DEany[DEany$Protein[!grepl("REV",DEany$Protein)],]#1547 observations


DEtable <- as.matrix(table(DEany$Protein))#the parenthetical should be a passed DE factor vector of networKin output
DEtable <- DEtable[DEtable!=0,,drop=F]#270 proteins


DEtable <- as.matrix(DEtable[row.names(DEtable) %in% row.names(BGtable),])#removing zeros and all factors not present in BG data
#note the use of the %in% statement to control for rownames in DeTable but not in background
#subset the background table in a similar way to ensure we are making the proper comparisons
BGtable <- as.matrix(BGtable[row.names(BGtable) %in% row.names(DEtable),])
NotDE <- BGtable-DEtable
facttemp <- as.factor(row.names(DEtable))###########not quite working as a factor
pvals <- c()


