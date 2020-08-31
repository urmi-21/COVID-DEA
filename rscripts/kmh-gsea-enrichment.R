# Performs various enrichment analyses of expression data
# from Black and White populations across TCGA and GTeX for
# Singh et. al.
#
# Kyle Hernandez, University of Chicago, kmhernan@uchicago.edu
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# Get entrez from symbol
to.entrez <- function(symbols) {
  bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)
}

# Load
df.gtex <- read.csv("/Users/kylehernandez/Projects/Other/COV-IRT/ancestry/eve-paper/gene-set/DEGs/gtex_de_results.csv")
df.tcga <- read.csv("/Users/kylehernandez/Projects/Other/COV-IRT/ancestry/eve-paper/gene-set/DEGs/tcga_de_results.csv")

# Out
odir <- "SupplementaryData/gsea_testout"

# Universe with entrez
uni.bitr <- to.entrez(df.tcga$Name)

# GSEA
## GTEX

### Set up inputs
gtex.entrez <- df.gtex %>%
  filter(Name %in% uni.bitr$SYMBOL) %>%
  left_join(uni.bitr, by=c("Name"="SYMBOL"))
gtex.fc <- gtex.entrez$logFC
names(gtex.fc) <- as.character(gtex.entrez$ENTREZID)
gtex.fc <- sort(gtex.fc, decreasing=TRUE)

### GSEA - GO Ontology

# Molecular Function
gtex.gse.mf <- gseGO(geneList=gtex.fc,
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",
                     by = "fgsea",
                     nPerm = 5000,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")
gtex.gse.mf.x <- setReadable(simplify(gtex.gse.mf), 'org.Hs.eg.db', 'ENTREZID')
gtex.gse.mf.cnet <- cnetplot(gtex.gse.mf.x, foldChange=gtex.fc, categorySize="pvalue", showCategory = 5)
# This kind of approach let's you manipulate the figure
gtex.gse.mf.cnet$scales$scales[[1]] <- scale_color_gradient(name="fold change", low="green", high="red", na.value = "yellow")
print(gtex.gse.mf.cnet)

# GSEA KEGG
gtex.gsea.kk <- gseKEGG(geneList=gtex.fc,
                        organism='hsa',
                        nPerm=5000,
                        pvalueCutoff = 0.05,
                        verbose = FALSE,
                        pAdjustMethod = "BH")
gtex.gsea.kk.x <- setReadable(gtex.gsea.kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(gtex.gsea.kk.x, foldChange=gtex.fc, categorySize="pvalue", showCategory = 5)
write.table(data.frame(gtex.gsea.kk.x), file=paste(odir, "gtex.pooled.gsea.kegg.tsv", sep="/"),
            sep="\t", row.names=FALSE, quote = FALSE)

## TCGA

tcga.entrez <- df.tcga %>%
  filter(Name %in% uni.bitr$SYMBOL) %>%
  left_join(uni.bitr, by=c("Name"="SYMBOL"))
tcga.fc <- tcga.entrez$logFC
names(tcga.fc) <- as.character(tcga.entrez$ENTREZID)
tcga.fc <- sort(tcga.fc, decreasing=TRUE)

# TCGA GSEA KEGG
tcga.gsea.kk <- gseKEGG(geneList=tcga.fc,
                        organism='hsa',
                        nPerm=5000,
                        by="fgsea",
                        pvalueCutoff = 0.05,
                        verbose = FALSE,
                        pAdjustMethod = "BH")
tcga.gsea.kk.x <- setReadable(tcga.gsea.kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(tcga.gsea.kk.x, foldChange=tcga.fc, categorySize="pvalue", showCategory = 5)
write.table(data.frame(tcga.gsea.kk.x), file=paste(odir, "tcga.pooled.gsea.kegg.tsv", sep="/"),
            sep="\t", row.names=FALSE, quote = FALSE)