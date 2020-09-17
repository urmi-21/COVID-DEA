library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(readr)



##########################################Functions###############################################
# Get entrez from symbol by KH
to.entrez <- function(symbols) {
  bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)
}

doKegg <- function(file,dir,outDir,suffix) {
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
  set.seed(21);
  gtex.gsea.kk <- gseKEGG(geneList=gtex.fc,
                          organism='hsa',
                          nPerm=5000,
                          pvalueCutoff = 0.5,
                          verbose = FALSE,
                          pAdjustMethod = "BH",
                          seed = T,
                          by = "fgsea"
                          )
  gtex.gsea.kk.x <- setReadable(gtex.gsea.kk, 'org.Hs.eg.db', 'ENTREZID')
  
  gtex.gsea.kk.cnet<-cnetplot(gtex.gsea.kk.x, foldChange=gtex.fc, categorySize="pvalue", showCategory = 5,
                              node_label=T)
  
  p1<-gtex.gsea.kk.cnet+ 
    scale_color_gradient2(name="logFC",low="grey",high="blue",midpoint = 0,na.value = "yellow")+
    guides( size = FALSE)+ 
    theme(legend.position = c(0.9, 0.2))+theme(text = element_text(size=15)) 
  #save plot 
  ggsave(paste(outDir,.Platform$file.sep,thisname,suffix,"_cnet.png",sep = ""),p1,width = 12,units = "in",dpi = 300)
  
  p2<-ridgeplot(gtex.gsea.kk)
  #save plot to file
  ggsave(paste(outDir,.Platform$file.sep,thisname,suffix,"_ridge.png",sep = ""),p2,width = 15,units = "in",dpi = 300)
  
  #write results
  write.table(data.frame(gtex.gsea.kk.x), file=paste(outDir,.Platform$file.sep,thisname,suffix,"_kegg.tsv",sep = ""),
              sep="\t", row.names=FALSE, quote = FALSE)
  
  #return
  return(data.frame(gtex.gsea.kk.x))
  
}


compare_topterms <- function(mwres,lmres,topn,name) {
  commonterms <- intersect(mwres$Description,lmres$Description)
  termsmw<-setdiff(mwres$Description,lmres$Description)
  termslm<-setdiff(lmres$Description,mwres$Description)
  
  terms_mw_top20 <- head(mwres$Description,topn)
  terms_lm_top20 <- head(lmres$Description,topn)
  commonterms_top20 <- intersect(terms_mw_top20,terms_lm_top20)
  termsmw_top20<-setdiff(terms_mw_top20,terms_lm_top20)
  termslm_top20<-setdiff(terms_lm_top20,terms_mw_top20)
  
  
  #write df
  results<-data.frame(Total_terms_MW=integer(),
                      Total_terms_Limma=integer(),
                      Intersection_All=integer(),
                      Only_MW=integer(),
                      Only_Limma=integer(),
                      Intersection_topn=integer(),
                      Only_MW_topn=integer(),
                      Only_Limma_topn=integer(),
                      Terms_Intersect_topn=character(),
                      Terms_Only_MW_topn=character(),
                      Terms_Only_Limma_topn=character(),
                      Terms_MW_topn=character(),
                      Terms_Limma_topn=character(),
                      Terms_Intersect_ALL=character(),
                      Terms_Only_MW=character(),
                      Terms_Only_Limma=character(),
                      stringsAsFactors=FALSE) 
  
  results[1,c('Total_terms_MW')]=c(length(mwres$Description))
  results[1,c('Total_terms_Limma')]=c(length(lmres$Description))
  results[1,c('Intersection_All')]=c(length(commonterms))
  results[1,c('Only_MW')]=c(length(termsmw))
  results[1,c('Only_Limma')]=c(length(termslm))
  results[1,c('Intersection_topn')]=c(length(commonterms_top20))
  results[1,c('Only_MW_topn')]=c(length(termsmw_top20))
  results[1,c('Only_Limma_topn')]=c(length(termslm_top20))
  
  row=1
  for(x in commonterms_top20){
    results[row,c('Terms_Intersect_topn')]<-x
    row=row+1
  }
  
  row=1
  for(x in termsmw_top20){
    results[row,c('Terms_Only_MW_topn')]<-x
    row=row+1
  }
  
  row=1
  for(x in termslm_top20){
    results[row,c('Terms_Only_Limma_topn')]<-x
    row=row+1
  }
  
  row=1
  for(x in terms_mw_top20){
    results[row,c('Terms_MW_topn')]<-x
    row=row+1
  }

  row=1
  for(x in terms_lm_top20){
    results[row,c('Terms_Limma_topn')]<-x
    row=row+1
  }
  
  row=1
  for(x in commonterms){
    results[row,c('Terms_Intersect_ALL')]<-x
    row=row+1
  }
  
  row=1
  for(x in termsmw){
    results[row,c('Terms_Only_MW')]<-x
    row=row+1
  }
  
  
  row=1
  for(x in termslm){
    results[row,c('Terms_Only_Limma')]<-x
    row=row+1
  }
  
  
  #remove NA
  results <- sapply(results, as.character)
  results[is.na(results)] <- ""
  
  #rename col headers
  #colnames(results)
  colnames(results)[6]<-paste0('Intersection_top_',topn)
  colnames(results)[7]<-paste0('Only_MW_top_',topn)
  colnames(results)[8]<-paste0('Only_Limma_top_',topn)
  colnames(results)[9]<-paste0('Terms_Intersect_top_',topn)
  colnames(results)[10]<-paste0('Terms_Only_MW_top_',topn)
  colnames(results)[11]<-paste0('Terms_Only_Limma_top_',topn)
  colnames(results)[12]<-paste0('Terms_MW_top_',topn)
  colnames(results)[13]<-paste0('Terms_Limma_top_',topn)
  
  
  
  #wrtite to file
  write.table(results, file=paste(outdir,.Platform$file.sep,name,"_compare.tsv",sep = ""),
              sep="\t", row.names=FALSE, quote = FALSE)
  
  }





#########################################END F############################################################


#tissuewise
limmadir<-"limma_tissuewise"
mwdir<-"MW_tissuewise"
outdir<-"mw_limma_compar"

infilesmw<-list.files(path = mwdir, pattern= "_mw.txt", recursive = TRUE)

infileslm<-list.files(path = limmadir, pattern= "_limmaOut.tsv", recursive = TRUE)

#dir="."
#infiles<-c("LSttest_BRCAsubtypes/BRCA_LStestOut.tsv")
#infiles<-c("limma_BRCAsubtypes/BRCA_limmaOut.tsv")
#doKegg(infiles[1],dir=dir,outDir = outdir)

smallmw<-head(infilesmw,3)
smalllm<-head(infileslm,3)
#make plots for MW
mwout<-apply(as.matrix(infilesmw),1,doKegg,dir=mwdir,outDir = outdir,suffix="_MW")
#make plots for limma
limmaout<-apply(as.matrix(infileslm),1,doKegg,dir=limmadir,outDir = outdir,suffix="_LM")




#mwres<-mwout[[1]]
#lmres<-limmaout[[1]]
#compare results
#compare_topterms <- function(mwres,lmres,topn,name) {
n=1
for (f in infilesmw){
  name<-strsplit(f,'_')[[1]][1]
  mwres<-mwout[[n]]
  lmres<-limmaout[[n]]
  compare_topterms(mwres,lmres,20,name)
  n=n+1
}






