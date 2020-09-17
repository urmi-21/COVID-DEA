library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(readr)


# Get entrez from symbol by KH
to.entrez <- function(symbols) {
  bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)
}

doKegg <- function(file,dir,outDir) {
  thisname<-strsplit(file,"_")[[1]][1]
  #read file
  df.gtex <- read_delim(paste(dir,.Platform$file.sep,file,sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
  # Universe with entrez
  uni.bitr <- to.entrez(df.gtex$Name)
  # GSEA
  ### Set up inputs
  gtex.entrez <- df.gtex %>%
    filter(Name %in% uni.bitr$SYMBOL) %>%
    left_join(uni.bitr, by=c("Name"="SYMBOL"))
  
  gtex.fc <- gtex.entrez$logFC
  names(gtex.fc) <- as.character(gtex.entrez$ENTREZID)
  gtex.fc <- sort(gtex.fc, decreasing=TRUE)
  # GSEA KEGG
  gtex.gsea.kk <- gseKEGG(geneList=gtex.fc,
                          organism='hsa',
                          nPerm=5000,
                          pvalueCutoff = 0.05,
                          verbose = FALSE,
                          pAdjustMethod = "BH")
  gtex.gsea.kk.x <- setReadable(gtex.gsea.kk, 'org.Hs.eg.db', 'ENTREZID')
  
  gtex.gsea.kk.cnet<-cnetplot(gtex.gsea.kk.x, foldChange=gtex.fc, categorySize="pvalue", showCategory = 5,
                              node_label=T)
  
  p1<-gtex.gsea.kk.cnet+ 
    scale_color_gradient2(name="logFC",low="grey",high="blue",midpoint = 0,na.value = "yellow")+
    guides( size = FALSE)+ 
    theme(legend.position = c(0.9, 0.2))+theme(text = element_text(size=15)) 
  #save plot 
  ggsave(paste(outDir,.Platform$file.sep,thisname,"_cnet.png",sep = ""),p1)
  
  p2<-ridgeplot(gtex.gsea.kk)
  #save plot to file
  ggsave(paste(outDir,.Platform$file.sep,thisname,"_ridge.png",sep = ""),p2)
  
  #write results
  write.table(data.frame(gtex.gsea.kk.x), file=paste(outDir,.Platform$file.sep,thisname,"_kegg.tsv",sep = ""),
              sep="\t", row.names=FALSE, quote = FALSE)
  
}


#tissuewise
dir<-"D:/urmi/docs/humandata/covid/bw_de"
outdir<-"tissuewise_kegg"
infiles<-list.files(path = dir, pattern= "_bw.txt", recursive = TRUE)

dir="."
#infiles<-c("LSttest_BRCAsubtypes/BRCA_LStestOut.tsv")
infiles<-c("limma_BRCAsubtypes/BRCA_limmaOut.tsv")
doKegg(infiles[1],dir=dir,outDir = outdir)

apply(as.matrix(infiles),1,doKegg,dir=dir,outDir = outdir)
