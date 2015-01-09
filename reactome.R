source("http://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
biocLite("ReactomePA")
require(DOSE)
##run through a demo
browseVignettes("ReactomePA")


# here is the demo
require(DOSE)
data(geneList)#from a cancer project listing ratios of expression across two conditions (using hypergeometric mean)
de <- names(geneList)[abs(geneList) > 1]
head(de)
require(ReactomePA)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)#see help file. uses hypergeometric distribution.
#note the universe background set and the cutoff used. I would want to use q values
head(summary(x))

barplot(x, showCategory=8)#showCategory is the number of categories to show

#output to data frame for further analysis
output <- as.data.frame(summary(x))

enrichMap(x, fixed =F)#not sure what this is. I think it links categories if they have shared genes.maybe useful

#shows the interconnections between categories and the fold change. I think reactome may be more useful in the end but perhaps interesting.How to export and play with this in cytoscape?
cnetplot(x, categorySize = "pvalue", foldChange = geneList)#unclear how to limit categories and play with output display
cnetplot(x, categorySize = "geneNum", foldChange = geneList)#unclear how to limit categories and play with output display
cnetplot(x, categorySize = "count", foldChange = geneList)#unclear how to limit categories and play with output display
cnetplot(x, showCategory = 1, categorySize = "pvalue", foldChange = geneList)#1 category at a time
cnetplot(x, showCategory = 1, categorySize = "pvalue", foldChange = geneList, fixed =F)#1 category at a time and open for layout editing!

#cluster pathway comparison visualization (can be used for different DE subset comparisons)#quite cool
require(clusterProfiler)
data(gcSample)
str(gcSample)#unclear what type of data this is
res <- compareCluster(gcSample, fun = "enrichPathway", universe = names(geneList))
res <- compareCluster(gcSample, fun = "enrichPathway")
plot(res)


##and now GSEA as opposed to pathway analysis.
y <- gsePathway(geneList, nPerm = 100, minGSSize = 120,
                pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
res <- summary(y)
head(res)

enrichMap(y)

gseaplot(y, geneSetID = "1280215")#look at where a given geneset falls in ranked list of genes

#Below a pathway viewer. I think a better choice may be found within the rectome website or perhaps using another software option like pathview.
viewPathway("E2F mediated regulation of DNA replication",
            readable = TRUE, foldChange = geneList)

# enrichment vs GSEA vs fishers exact test with mtesting correction vs directionality
#fisher's exact test and hypergeotric test really the same 

#difference between GSEA and enrichment analysis (hypergeometric dist/Fisher's exact test style approaches) is that the GSEA uses all genes and tries to identify DE gene sets while the enrichment analysis looks at DE sets of genes to look for enrichment of GO/pathway or any other categorical-type of enrichment. Here they just use the fold change to calculate the results for GSEA.
