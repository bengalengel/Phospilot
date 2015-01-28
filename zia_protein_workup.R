##zia protein workup. goal is to get protein measurements for all 60 human samples. Normalize them against each other (median then quantile normalized) and then subtract (log scale) the relative protein concentration measurements from the relative phospho measurements. Issues will be in isoform quantification and batch correction across the protein measurements.

#these measurements will then be subjected to the same workflow. 

rm(list=ls(all=TRUE)) #start with empty workspace


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)
source("loadMQ.R")
source("ExpandPhos.R")
source("counts.R")
source("breakdown.R")

# load protein files with particular variables populated using "loadMQ"
protein <- load.MQ(directory = "C:/Users/Brett/Documents/Pilot/10_9_14/txt/", type = "protein")


# load protein files with particular variables populated using "loadMQ" at home
protein <- load.MQ(directory = "E:/My Documents/Pilot/10_9_14/txt/", type = "protein")


# remove contaminants and reverse database hits
phospho <- phospho[(phospho$Potential.contaminant != "+" & phospho$Reverse != "+"),]
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]

#expand for observation (multiplicity based analysis). All observations quantified in >=1 experiment.
multExpanded <- ExpandPhos(phospho)

# subset phospho to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75),]

# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
# Note there are some high probability proteins here given some proteins are idd by 20+ phosphopeptides.
# eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"),]

# Class 1 sites with each source of quantification for that site (singly/doubly/3+) explicitly accounted for 
multExpanded1 <- ExpandPhos(phospho1)

#make tables of basic counts of proteins and phosphopeptides (may have to update when performing the normalization)
phoscount(phospho,phospho1,multExpanded,multExpanded1)
proteincount(protein)

# make breakdown charts of phospho and protein overlap. This function produces barplots of number of sites/obs and proteins per experiment and
# cumulative over experiments. It also has extensive phospo info, including number of phospho per protein, multiplicity breakdown, and venns
breakdown(protein, phospho, multExpanded, cls=F)
breakdown(protein, phospho1, multExpanded1)
