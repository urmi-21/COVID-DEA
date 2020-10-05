library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(readr)

#Perform tissue/tumor-wise GSEA using MW and Limma DE genes and do comparison



##########################################Functions###############################################
# Get entrez from symbol by KH
to.entrez <- function(symbols) {
  bitr(symbols, "SYMBOL", "ENTREZID", org.Hs.eg.db)
}

doKegg <- function(file,dir,outDir,suffix) {
  thisname<-strsplit(file,"_")[[1]][1]
  print(paste0(thisname,"_",suffix))

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
                          pvalueCutoff = 0.1,
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
  print(paste(outDir,.Platform$file.sep,thisname,suffix,"_kegg.tsv",sep = ""))
  write.table(data.frame(gtex.gsea.kk.x), file=paste(outDir,.Platform$file.sep,thisname,suffix,"_kegg.tsv",sep = ""),
              sep="\t", row.names=FALSE, quote = FALSE)
  
  #return
  return(data.frame(gtex.gsea.kk.x))
  
}


compare_topterms <- function(mwres,lmres,topn,name) {
  
  #sort by NES
  #mwres <- mwres %>% arrange(desc(NES))
  #lmres <- lmres %>% arrange(desc(NES))
  
  #sort by abs value
  #mwres <- mwres %>% arrange(desc(abs(NES)))
  #lmres <- lmres %>% arrange(desc(abs(NES)))
  
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


compare_degenes <- function(name,mwdir,limmadir,outDir) {
 #compare tissue wise DE genes
 df.mw <- read_delim(paste(mwdir,.Platform$file.sep,paste0(name,"_mw.txt"),sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
 df.limma <- read_delim(paste(limmadir,.Platform$file.sep,paste0(name,"_limmaOut.tsv"),sep=""), "\t", escape_double = FALSE, trim_ws = TRUE)
 #sort by fc
 df.mw <- df.mw %>% arrange(desc(logFC))
 df.limma <- df.limma %>% arrange(desc(logFC))
 

 upmw<-df.mw %>% filter(`Adj pval` <0.05) %>% filter(logFC >= 1)
 dwnmw<-df.mw %>% filter(`Adj pval` <0.05) %>% filter(logFC <= -1)
 
 uplm<-df.limma %>% filter(`Adj pval` <0.05) %>% filter(logFC >= 0.5)
 dwnlm<-df.limma %>% filter(`Adj pval` <0.05) %>% filter(logFC <= -0.5)
 
 #top 50 genes
 N=100
 topupmw<-head(upmw,N)$Name
 topdwnmw<-tail(dwnmw,N)$Name
 
 topuplm<-head(uplm,N)$Name
 topdwnlm<-tail(dwnlm,N)$Name
 
 upintersect<-intersect(topupmw,topuplm)
 uponlymw<-setdiff(topupmw,topuplm)
 uponlylm<-setdiff(topuplm,topupmw)
 
 dwnintersect<-intersect(topdwnmw,topdwnlm)
 dwnonlymw<-setdiff(topdwnmw,topdwnlm)
 dwnonlylm<-setdiff(topdwnlm,topdwnmw)

 #write to file
 results<-data.frame(Up_MW=character(),
                     Up_limma=character(),
                     Up_Intersection_All=character(),
                     Up_Only_MW=integer(),
                     Up_Only_Limma=integer(),
                     Dwn_MW=character(),
                     Dwn_limma=character(),
                     Dwn_Intersection_All=character(),
                     Dwn_Only_MW=integer(),
                     Dwn_Only_Limma=integer(),
                     stringsAsFactors=FALSE) 
 
 row=1
 for(x in topupmw){
   results[row,c('Up_MW')]<-x
   row=row+1
 }
 
 row=1
 for(x in topuplm){
   results[row,c('Up_limma')]<-x
   row=row+1
 }
 
 row=1
 for(x in upintersect){
   results[row,c('Up_Intersection_All')]<-x
   row=row+1
 }
 
 row=1
 for(x in uponlymw){
   results[row,c('Up_Only_MW')]<-x
   row=row+1
 }

 row=1
 for(x in uponlylm){
   results[row,c('Up_Only_Limma')]<-x
   row=row+1
 }
 
 row=1
 for(x in topdwnmw){
   results[row,c('Dwn_MW')]<-x
   row=row+1
 }
 
 row=1
 for(x in topdwnlm){
   results[row,c('Dwn_Limma')]<-x
   row=row+1
 }
 
 row=1
 for(x in dwnintersect){
   results[row,c('Dwn_Intersection_All')]<-x
   row=row+1
 }
 
 row=1
 for(x in dwnonlymw){
   results[row,c('Dwn_Only_MW')]<-x
   row=row+1
 }
 
 row=1
 for(x in dwnonlylm){
   results[row,c('Dwn_Only_Limma')]<-x
   row=row+1
 }
 
 #remove NA
 results <- sapply(results, as.character)
 results[is.na(results)] <- ""
 
 #wrtite to file
 write.table(results, file=paste(outDir,.Platform$file.sep,name,"_DE_compare.tsv",sep = ""),
             sep="\t", row.names=FALSE, quote = FALSE)
   
}





#########################################END F############################################################


#dir with DE results
limmadir<-"limma_tissuewise"
#limmadir<-"limma_tissuewise_bmi"
mwdir<-"MW_tissuewise"
outdir<-"mw_limma_compare5"
#create outdir if not present
dir.create(file.path(outdir))

infilesmw<-list.files(path = mwdir, pattern= "_mw.txt", recursive = TRUE)

infileslm<-list.files(path = limmadir, pattern= "_limmaOut.tsv", recursive = TRUE)

#dir="."
#infiles<-c("LSttest_BRCAsubtypes/BRCA_LStestOut.tsv")
#infiles<-c("limma_BRCAsubtypes/BRCA_limmaOut.tsv")
#doKegg(infiles[1],dir=dir,outDir = outdir)
#smallmw<-head(infilesmw,3)
#smalllm<-head(infileslm,3)

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


names<-apply(as.matrix(infilesmw),1,function(file){return(strsplit(file,"_")[[1]][1])})
#name<-names[1]
#mwfiles<-infilesmw[1]
#lmfiles<-infileslm[1]

apply(as.matrix(names),1,compare_degenes,mwdir=mwdir,limmadir=limmadir,outDir = outdir)





