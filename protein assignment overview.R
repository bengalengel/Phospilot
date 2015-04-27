##There are 41 mappings (out of 18238) that are different between two approaches. 

# The first/original approach was the following.USING DE-ISOFORMED DATA.
# 1) Use the 'protein' column in the posphosite table to map to ANY of the proteins contained within the 'majority protein id' column from the proteingroups file derived from the proteome data.
# 2) If multiple protein groups are found, assign the protein group with the most ids (unique plus razor) to the phosphopeptide.

# The latest approach 
# 1) Find the phosphopeptide sequence in the proteome (using the FASTA db file)
# 2) Search the protein workup for these peptides
# 3) If multiple groups are found, assign the protein group with the most ids
# 4) Remove protein ids and majority protein ids if they don't contain the peptide of interest

#Weakness of approach 1

# 1) If a phosphopeptide is unique to an isoform this would result in incorrect mapping.
# 2) The majority protein ids reported may not contain the phosphopeptide sequence

#weaknesses of both approaches
# 1) The quantifications assigned to a phosphopeptide for normalization could arise from any of the proteins assigned to the group
# 2) Protein group assignment accuracy is limited by the depth of the data analysis (number of proteotypic peptides identified) and is IMPERFECT


# A spot check of the differences between the two approaches shows that the discrepancies are entirely a result of choosing the wrong isoforms for protein assignment

#approach number two is used because it doesn't suffer from isoform issues

# Issue 1 for common weaknesses is being resolved by removing proteins and majority proteins that don't contain the phosphopeptide and therefore should not be used for ANNOTATION PURPOSES. Since there is no way disambiguiate the contributions of protein group members to the protein group quantifications there is no way to adjust the quantifications. This should have no effect on protein FDR.