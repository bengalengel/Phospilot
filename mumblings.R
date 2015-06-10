
#add a +/- based on presense of ENSP ID within any of the "leading proteins" for confounded data and within "ppMajorityProteinIDs" for Zia workup normalized protein ids. Also pass 1) which of the proteins within the groups match for each phosphopeptide 2) homo or heterozygote and 3) snp number

# I am going to need a new dataframe for this...

#I can use mapply for the first part. but I want to merge all of the data from the SNPenrich DF together with the me df
1 For each cell line x, (x = 1:4) does any protein group member assigned to peptide z (z = 1:nrow(multexpanded1))
match the protein associated with snp y (y =1:nrow(SNPeffFinal))? (+/-)
2 which proteins within the group matched a snp?
3 For those proteins within the group that matched a snp list the snps that matched it... fuckedy duck

...OK 1st thing is see if there is an enrichment.

#things to note for the future are the positions of the phosphosite within the protein relative to the SNP, domain position within the protein, proximity of a snp to a phosphorylated residue, proximity of the snp toward a residue that has been annotated as being phosphorylated. These things may be difficult to handle due to annotation issues. 

#feels too complicated for ifelse, I will use mapply with a named function 'MatchProtSNPeff' I don't want to iterate over 'SampleSNP'
MatchProtSNPeff <- function(queryproteins, SampleSNP){
  #queryproteins are proteins assigned to the peptide and sampleSNP is the 'genotype' for all the snps for that sample. Note that there can be 1 to many relationships here...
  if(any(unlist(strsplit(as.character(queryproteins), ";")) %in% SNPeffFinal$peptide) &&  any(hapTypes %in% SampleSNP)){
    "+"
  }else{
    "-"
  }
  
  
  #I can't do this with variable length vectors. its not clear I should do this anyway. MANY TO MANY RELATIONSHIP STRUCTURE HERE!!
  #I want 
  #For enrichment analysis in diffphos
  1) does this peptide match any snp for any of the samples (BY WAY OF PROTEIN)?
2) How many snps match this peptide? (to show distribution of snps/protein group in my samples. Bias in enricment should be accounted for by the background)
#For effect size comparisons between lines
3) For a given snp, which lines have it?
4) For each line that has the snp, what is its genotype?

# I am also interested in specific snp properties such as:
#for each protein/isoform within the group assigned to a peptide

1) Does it result in the creation/deletion of phos/ser/thr?
2) Is it adjacent to the measured phosphopeptide?
3) Is it equivalent to the measured phosphopeptide?
4) Is it within a domain?

#The above comes with the caveat that positional effects for the snp may manifest differently depending on the isoform. If isoforms are ambiguous a flag should be passed to indicate this.
