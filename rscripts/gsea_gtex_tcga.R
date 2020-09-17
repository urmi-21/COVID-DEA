# Performs various enrichment analyses of expression data
# from Black and White populations across TCGA and GTeX for
# Singh et. al.
#
# Kyle Hernandez, University of Chicago, kmhernan@uchicago.edu
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
# Get entrez from symbol
to.entrez <- function(symbols) {
  bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)
}

# Load
df.gtex <- read_delim("gtex_de.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
df.tcga <- read_delim("tcga_de.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
# Out
odir <- "gtex_tcga_gsea"

# Universe with entrez
uni.bitr <- to.entrez(df.tcga$Name)

### Set up inputs
gtex.entrez <- df.gtex %>%
  filter(Name %in% uni.bitr$SYMBOL) %>%
  left_join(uni.bitr, by=c("Name"="SYMBOL"))

gtex.fc <- gtex.entrez$logFC
names(gtex.fc) <- as.character(gtex.entrez$ENTREZID)
gtex.fc <- sort(gtex.fc, decreasing=TRUE)


# GSEA KEGG
set.seed(22);
gtex.gsea.kk <- gseKEGG(geneList=gtex.fc,
                        organism='hsa',
                        nPerm=5000,
                        pvalueCutoff = 0.05,
                        verbose = FALSE,
                        pAdjustMethod = "BH")
gtex.gsea.kk.x <- setReadable(gtex.gsea.kk, 'org.Hs.eg.db', 'ENTREZID')

gtex.gsea.kk.cnet<-cnetplot(gtex.gsea.kk.x, foldChange=gtex.fc, categorySize="pvalue", showCategory = 5,
                            node_label=T)

gtex.gsea.kk.cnet+ 
  scale_color_gradient2(name="logFC",low="grey",high="blue",midpoint = -1.,na.value = "yellow")+
  guides( size = FALSE)+ 
  theme(legend.position = c(0.9, 0.2))+theme(text = element_text(size=15)) 


ridgeplot(gtex.gsea.kk)

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
tcga.gsea.kk.cnet<-cnetplot(tcga.gsea.kk.x, foldChange=tcga.fc, categorySize="pvalue", showCategory = 5)

tcga.gsea.kk.cnet+ 
  scale_color_gradient2(name="logFC",low="red",high="grey",midpoint = 0,na.value = "yellow")+
  guides( size = FALSE)+ 
  theme(legend.position = c(0.9, 0.1))+theme(text = element_text(size=15)) 

write.table(data.frame(tcga.gsea.kk.x), file=paste(odir, "tcga.pooled.gsea.kegg.tsv", sep="/"),
            sep="\t", row.names=FALSE, quote = FALSE)






