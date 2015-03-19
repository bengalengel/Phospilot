# adding the annotations and the pudding, with the jellow (cosby voice)...too soon?

#assign to multExpanded1 (proteins - that is the majority proteins from each of the protein groups assigned to a phosphopeptide. Perhaps the number of isoforms has to do with the lack of tie-breaker approach used at the protein group level. That is, simply the highest ranked protein is returned.) GO terms using R annotation packages


## ----echo=FALSE----------------------------------------------------------
suppressPackageStartupMessages(library(org.Hs.eg.db))

#load the human gene centric OrgDb package.
## ------------------------------------------------------------------------
library(org.Hs.eg.db)

#I can create my own OrganismDb package or use an existing one with diverse annotations

#annotationhub packages allow me to get access to a LARGE variety of external databases.

library(AnnotationHub)

#I have to create my own ah object
ah = AnnotationHub()#must be online. This is pretty powerful


#biomaRt exposes a huge family of online annotation resources called marts. Load the backage and decide which 'mart' I want to use. useMart() selects the annoation mart.

library("biomaRt")
listMarts()#need the interwebs for this. Data is current and should work from this locally once it is downloaded so that things are reproducible.

#anyway after I choose a mart I need to decide on the dataset. use the listDatasets(martobject) function
#use the useMart() function again to decide on a dataset.

#using a combination of filters and attributes I can retrieve stuff using 'getBM'

#I NEED TO PAY ATTENTION TO THE GENOME BUILD WHEN RETRIEVING ANNOTATIONS!!

#BSgenome objects exist for extracting sequence information
library(AnnotationDbi)
help(package="AnnotationDbi")

!!!!!!!!!!!!!!!!!!!!!!
#variantannotation package to assess if something is within a coding region or what. I should produce a figure using this package. 

# OrganismDb packages are named for the species they represent (such as the Homo.sapiens package). These packages contain references to other key annotations packages and can thus represent all the underlying data as if it were coming from one place. So for example, the Homo.sapiens package can allow you to retrieve data about the ranges of a genes transcripts at the same time that you extract it's gene name because it represents both a the transcriptome and the relevant org package for Homo sapiens. These can be generated using functions in the OrganismDbi package if you have specific packages that you want to link together.












