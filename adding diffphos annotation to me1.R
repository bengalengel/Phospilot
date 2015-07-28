
#add annotations for confounded data with blocking and technical duplicate consensus correlation information

#add annotation to multexpanded DF
multExpanded1$SubtoDEConfounded = ifelse(multExpanded1$idmult %in% row.names(adata),"+","-")

#add F test values to the table
multExpanded1$globalFsigConfounded = ifelse(multExpanded1$idmult %in% row.names(sigFvalsCombat),"+","-")


#phosprep as a covariate

#add annotation to multexpanded DF
multExpanded1$SubtoDEPhosProt = ifelse(multExpanded1$idmult %in% row.names(PhosProt),"+","-")

#add F test values to the table
multExpanded1$globalFsigPhosProt = ifelse(multExpanded1$idmult %in% row.names(sigFvalsPhosPrepProt),"+","-")


#GelPrep as a covariate

#add annotation to multexpanded DF
multExpanded1$SubtoDEGelProt = ifelse(multExpanded1$idmult %in% row.names(PhosProt2),"+","-")

#add F test values to the table
multExpanded1$globalFsigGelProt = ifelse(multExpanded1$idmult %in% row.names(sigFvalsGelPrepProt),"+","-")






