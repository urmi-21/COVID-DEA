library(readr)
library(dplyr)
library(diptest)
library(plyr)
library(scales)
library(data.table)
library(ggplot2)
library(ggthemes)
library(ggpubr)


#This script takes 3 arguments
#arg1: Path to input file (this file is created by MOG)
#arg2: Path to the metadatafile (used in the MOG project)
#arg3: Path to the output directory (specified by the MOG user in MOG)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Invalid Arguments", call.=FALSE)
} 

#read data
infile<-args[1]
metadataFile<-args[2]
outDir<-args[3]


#read metadata
metadata<- read_delim(metadataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(metadata)[which(colnames(metadata)=="portions.analytes.aliquots.submitter_id")]='SID'
#filter metadata
races<-c("black or african american","white")
diseases<-c('BRCA','COAD','KIRC','KIRP','LUAD','LUSC','THCA','UCEC')
metadata<-metadata%>%filter(clinical.race %in% races) %>% filter(clinical.disease %in% diseases)

#read data
data<-read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)

allGenes<-data$Name
#logFile = paste(outDir,.Platform$file.sep,'log_file.txt',sep = "")



#cat("Start log ", file=logFile, append=FALSE, sep = "\n")

#perform dip test for all genes
dipTestResults<-data.frame(Disease=character(),
                           Sample=character(), 
                           Gene=character(),
                           Data=character(),
                           Size=integer(),
                           Method=character(),
                           Statistic=double(),
                           pval=double(),
                           stringsAsFactors=FALSE) 

#cat("Start DT ... ", file=logFile, append=TRUE, sep = "\n")


for(thisDis in diseases){
  #print(thisDis)
  
  for (geneind in 1:length(allGenes)) {
    thisGene<-allGenes[geneind]
    #print(thisGene)
    for (thisRace in races) {
      #print(thisRace)
      thisSamps<-metadata %>% filter(clinical.race == thisRace) %>% filter(clinical.disease == thisDis)
      thisData<-data[geneind,thisSamps$SID]    
      tmp<-dip.test(log2(as.numeric(thisData)+1))
      dipTestResults[nrow(dipTestResults) + 1,] = c(thisDis,thisRace,thisGene,tmp$data.name,tmp$nobs,tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))

    }
    
  }
}



#cat("Finish DT ... ", file=logFile, append=TRUE, sep = "\n")

#perform KS test
#perform dip test for all genes
KSTestResults<-data.frame(Disease=character(),
                          Gene=character(),
                          Data=character(),
                          Size_grp1=integer(),
                          Size_grp2=integer(),
                          Method=character(),
                          Statistic=double(),
                          pval=double(),
                          stringsAsFactors=FALSE) 

for(thisDis in diseases){
  #print(thisDis)
  for (geneind in 1:length(allGenes)) {
    thisGene=allGenes[geneind]
    #print(thisGene)
    
    blackSamps<-metadata %>% filter(clinical.race == races[1]) %>% filter(clinical.disease == thisDis)
    whiteSamps<-metadata %>% filter(clinical.race == races[2]) %>% filter(clinical.disease == thisDis)
    blackData<-data[geneind,blackSamps$SID]
    whiteData<-data[geneind,whiteSamps$SID]
    #ks.test
    tmp<-ks.test(log2(as.numeric(blackData)+1),log2(as.numeric(whiteData)+1))
    KSTestResults[nrow(KSTestResults) + 1,] = c(thisDis,thisGene,tmp$data.name,length(blackData),length(whiteData),tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
    
  }
}

#cat("Finish KS ... ", file=logFile, append=TRUE, sep = "\n")

#compute foldchange
FCResults<-data.frame(Disease=character(),
                      Gene=character(),
                      Size_grp1=integer(),
                      Size_grp2=integer(),
                      Mean_grp1=integer(),
                      Mean_grp2=integer(),
                      LogFC=double(),
                      FC=double(),
                      stringsAsFactors=FALSE) 

for(thisDis in diseases){
  #print(thisDis)
  for (geneind in 1:length(allGenes)) {
    thisGene=allGenes[geneind]
    #print(thisGene)
    
    blackSamps<-metadata %>% filter(clinical.race == races[1]) %>% filter(clinical.disease == thisDis)
    whiteSamps<-metadata %>% filter(clinical.race == races[2]) %>% filter(clinical.disease == thisDis)
    blackData<-data[geneind,blackSamps$SID]
    whiteData<-data[geneind,whiteSamps$SID]
    
    #compute log FC
    logfc=mean(log2(as.numeric(blackData)+1))-mean(log2(as.numeric(whiteData)+1))
    fc=2^logfc
    m1=mean(as.numeric(blackData)+1)
    m2=mean(as.numeric(whiteData)+1)
    
    FCResults[nrow(FCResults) + 1,] = c(thisDis,thisGene,length(blackData),length(whiteData),m1,m2,
                                        logfc,fc)
  }
}

#cat("Finish FC ... ", file=logFile, append=TRUE, sep = "\n")

#save results
#write_tsv(dipTestResults,paste(outDir,.Platform$file.sep,'diptest.tsv',sep = ""))
#write_tsv(KSTestResults,paste(outDir,.Platform$file.sep,'KStest.tsv',sep = ""))

#cat("Finish save results... ", file=logFile, append=TRUE, sep = "\n")

#Read file containing significant MWtest
#MOG_MW <- read_delim("MOG_MW.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
MOG_MW<-data.frame(Disease=character(),
                   Gene=character(),
                   FDR=double(),
                   stringsAsFactors=FALSE)
#add data manually, right now arguments are not flexible
#only add significant FDR values. Full data is provided as SI
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','CCL3',0.004738966)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','CCL3L3',1.75E-12)
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','CCL3L3',0.015382274)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRP','CCL3L3',0.011288233)
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','CXCL9',0.00521798)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','CXCL9',0.007493509)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRP','CXCL9',0.031575819)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','CXCL10',0.009641767)
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','CXCL10',6.51E-04)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','CXCL10',0.029443046)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','F8A1',0)
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','F8A1',3.05E-06)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','F8A1',1.98E-07)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRP','F8A1',7.11E-09)
MOG_MW[nrow(MOG_MW) + 1,] = c('LUAD','F8A1',6.08E-06)
MOG_MW[nrow(MOG_MW) + 1,] = c('LUSC','F8A1',0.042371112)
MOG_MW[nrow(MOG_MW) + 1,] = c('THCA','F8A1',1.17E-04)
MOG_MW[nrow(MOG_MW) + 1,] = c('UCEC','F8A1',3.86E-06)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','F8A2',0)
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','F8A2',1.86E-05)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','F8A2',1.63E-08)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRP','F8A2',1.12E-08)
MOG_MW[nrow(MOG_MW) + 1,] = c('LUAD','F8A2',5.75E-07)
MOG_MW[nrow(MOG_MW) + 1,] = c('LUSC','F8A2',0.009591512)
MOG_MW[nrow(MOG_MW) + 1,] = c('THCA','F8A2',7.94E-04)
MOG_MW[nrow(MOG_MW) + 1,] = c('UCEC','F8A2',1.99E-04)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','F8A3',0)
MOG_MW[nrow(MOG_MW) + 1,] = c('COAD','F8A3',0.001375925)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','F8A3',1.25E-07)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRP','F8A3',1.04E-07)
MOG_MW[nrow(MOG_MW) + 1,] = c('LUAD','F8A3',3.55E-08)
MOG_MW[nrow(MOG_MW) + 1,] = c('LUSC','F8A3',6.60E-04)
MOG_MW[nrow(MOG_MW) + 1,] = c('UCEC','F8A3',0.006165061)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','GSTM1',4.11E-05)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','GSTM1',0.002426698)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','IL6ST',1.63E-13)
MOG_MW[nrow(MOG_MW) + 1,] = c('KIRC','IL6ST',0.00301579)
MOG_MW[nrow(MOG_MW) + 1,] = c('BRCA','TMPRSS6',1.23E-04)
MOG_MW$FDR<-as.numeric(MOG_MW$FDR)

#cat("Finish loading MW results... ", file=logFile, append=TRUE, sep = "\n")
                              
#create plots
#do for each row in data
makeplot <- function(x) {
  
  #return(as.data.frame(x,col.names =NULL,stringsAsFactors=FALSE))
  tdata<-as.data.frame(x,col.names =NULL,stringsAsFactors=FALSE)
  #return(tdata)
  tdata$SID<-rownames(tdata)
  names(tdata)<-tdata[1,]
  names(tdata)[2]<-'SID'
  tdata <- tdata[-1, ]
  rownames(tdata)<-NULL
  #return(tdata)
  genename<-colnames(tdata)[1]
  #print(paste('dim:',dim(tdata),' gene:',genename))
  
  #join with metadata
  tdata<-merge(tdata,metadata,all.x = TRUE,by='SID')
  #return(tdata)
  testdata2<- tdata %>% select('clinical.disease',genename,'clinical.race')
  colnames(testdata2)[1]='Disease'
  colnames(testdata2)[3]='Race'
  
  #print(paste('dim:',dim(testdata2)))
  
  #subset race
  testdata2<-testdata2 %>% filter(Race %in% c('white','black or african american'))
  testdata2$Race[which(testdata2$Race=='white')]<-'EA'
  testdata2$Race[which(testdata2$Race=='black or african american')]<-'AA'
  testdata2$Race<-as.factor(testdata2$Race)
  testdata2[,which(colnames(testdata2)==genename)]=as.numeric(testdata2[,which(colnames(testdata2)==genename)]) 
  #return(testdata2)
  
  #filter
  #testdata2<-testdata2%>%filter(Disease %in% c('COAD','BRCA','KIRC','KIRP','LUAD','LUSC','THCA','UCEC'))
  
  medians<-aggregate( get(genename) ~ Disease+Race, testdata2, median )
  colnames(medians)[3]<-'median'
  means<-aggregate( get(genename) ~ Disease+Race, testdata2, mean )
  colnames(means)[3]<-'mean'
  #trick to get legend
  emptylines<-means %>% mutate(Race=ifelse(Race=='AA','HD','KS'))
  #create df for three classes KE, HD and MW
  emptylines<-rbind(emptylines,means[1:8,]%>%mutate(Race='MW'))
  emptylines<-emptylines %>% mutate(mean=0) 
  colnames(emptylines)[which(colnames(emptylines)=='Race')]<-'Test'
  
  #plot violin plot
  dodge <- position_dodge(width = 0.9)
  
  ###Do calculations to add significance
  testdata2$Race<- as.factor(testdata2$Race)
  testdata2$Disease<-as.factor(testdata2$Disease)
  
  tdst<-compare_means(
    as.formula(paste(genename,"~ Race")), data = testdata2, group.by = "Disease",
    method = "t.test", ref.group = "AA"
  )
  
  #print(paste("Done",genename,sep="::"))
  #get position to display significance
  
  #maxVal<-max(testdata2[,genename])
  maxVals<-aggregate(as.formula(paste(genename,"~ Disease")), data = testdata2, max)
  colnames(maxVals)[2]<-"y.position"
  #increment_pct<-0.5
  #increment<-0
  #dip <- join(tdst,maxVals)
  #KS<-dip%>% mutate(y.position=y.position+increment_pct*y.position)
  #MW<-dip%>% mutate(y.position=y.position+(2*increment_pct*y.position))
  #FCval<-dip%>% mutate(y.position=y.position+(3*increment_pct*y.position))
  dip <- join(tdst,maxVals)
  KS<-dip
  MW<-dip
  FCval<-dip
  #order to add: diptest, KS, MW, FC
  
  #add diptest res
  diptmp<-dipTestResults%>%filter(Gene == genename)%>%filter(pval<0.05)
  dip<-dip%>%mutate(sigAA=" ")%>%mutate(sigEA=" ")  
  for(i in rownames(diptmp)){
    thisD<-diptmp[i,"Disease"]
    thisR<-diptmp[i,"Sample"]
    if (thisR=="black or african american"){
      dip[which(dip$Disease == thisD),"sigAA"]<-"*"
    }
    if (thisR=="white"){
      dip[which(dip$Disease == thisD),"sigEA"]<-"*"
    }
    
  }
  #join two sig in one column for display
  dip <- dip%>%mutate(p.signif=paste(sigAA,sigEA,sep="    "))
  
  #Add KS
  kstmp<-KSTestResults%>%filter(Gene == genename)%>%filter(pval<0.05)
  kssig<-unique(kstmp$Disease)
  KS<-KS%>%mutate(p.signif= ifelse(Disease %in% kssig , "*", ""))
  
  #add MW
  mwtmp<-MOG_MW%>%filter(Gene == genename)%>%filter(FDR<0.05)
  mwsig<-unique(mwtmp$Disease)
  MW<-MW%>%mutate(p.signif= ifelse(Disease %in% mwsig , "*", ""))
  
  #addFC
  fctmp<-FCResults%>%filter(Gene==genename)
  FCval<-join(FCval,fctmp)
  FCval$FC<-format(round(as.numeric(FCval$FC), 2), nsmall = 2)
  FCval<-FCval%>% mutate(FC=paste("FC:",FC,sep=" "))
  
  #handle position
  KS$y.position<-dip$y.position
  MW$y.position<-dip$y.position
  FCval$y.position<-dip$y.position
  
  increment_pct<-0.5
  #manage positions for pv labels
  for(i in rownames(dip)){
    thisD<-dip[i,"Disease"]
    #print(thisD)
    signOrder<-c(dip[i,"p.signif"],KS[i,"p.signif"],MW[i,"p.signif"],FCval[i,"p.signif"])
    #print(signOrder)
    thisStartVal<-dip[i,"y.position"]
    finalPos<-c(thisStartVal,thisStartVal,thisStartVal,thisStartVal)
    
    for(s in 1:length(signOrder)){
      thisSig<-signOrder[s]
      if(grepl("*", thisSig, fixed = TRUE)){
        for(k in s+1:length(signOrder)){
          finalPos[k]<-finalPos[k]+increment_pct*finalPos[k]
        }
      }
    }
    
    #update positions
    dip[i,"y.position"]<-finalPos[1]
    KS[i,"y.position"]<-finalPos[2]
    MW[i,"y.position"]<-finalPos[3]
    FCval[i,"y.position"]<-finalPos[4]
    
  }
  FCval<-FCval%>%mutate(y.position=y.position+0.5*y.position)
  
  #plot
  #title
  #ann_text<-data.frame(Type = c('a','b','c','d','e','f','g','h'), lbl = c('', '', '','','','','','sdas'))
  
  ggplot(testdata2, aes_string(x='Disease',  y=genename))+
    #annotate("text",  x=Inf, y = Inf, label = "Some text", vjust=1, hjust=1,colour = "red", size = 1.5)+
    geom_violin(aes(fill=Race),position = dodge,alpha = 0.1,show.legend = FALSE)+
    geom_boxplot(aes(fill=Race),width=.1, outlier.colour=NA, position = dodge,show.legend = FALSE)+
    stat_pvalue_manual(dip,x='Disease',color = "Blue",label = "p.signif",size = 5)+
    stat_pvalue_manual(KS,x='Disease',color = "firebrick3",label = "p.signif",size = 5)+
    stat_pvalue_manual(MW,x='Disease',color = "seagreen",label = "p.signif",size = 5)+
    stat_pvalue_manual(FCval,x='Disease',color = "Black",label = "FC",size = 2.5)+
    scale_y_log10()+
    theme(legend.key=element_blank(),legend.position="top",
          legend.title=element_blank(),legend.text = element_text(size=7),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(size = 8,face = "bold"),
          axis.text.y = element_text(size = 6,face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), axis.title=element_text(size=7,face="bold"))+
    
    ylab(paste(genename,"(FPKM)"))+xlab("")+
    geom_hline(data = means, aes(yintercept = mean, color = Race))+ 
    guides(colour=FALSE)+guides(color = guide_legend(override.aes = list(size = 3)))+
    #geom_hline(data = emptylines,alpha=0, aes(yintercept = mean, color = Test))+ 
    #scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
    #scale_colour_manual(name = 'the colour')+
    #guides(colour=FALSE)+guides(color = guide_legend(override.aes = list(size = 1)))+
    #geom_hline(data = medians, aes(yintercept = median))+
    #geom_text(data = ann_text, aes(x = Inf, y = Inf, label = lbl), hjust = 1,vjust=1)+
    facet_wrap( ~ Disease, scales="free",ncol = 8)+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),panel.spacing = unit(0.1, "lines"),panel.margin.x=unit(0.1, "lines")
    )#+annotate("text",  x=Inf, y = Inf, label = "Some text", vjust=1, hjust=1,colour = "red", size = 1.5)
  
  
  #print(paste('saving:',genename))
  #plot
  ggsave(paste(outDir,.Platform$file.sep,genename,".png",sep = ""),width = 8, height = 3, units = "in",dpi = 300)
  
}

#cat("Finish plot function... ", file=logFile, append=TRUE, sep = "\n")

#td<-apply(data[1,], 1, makeplot)
apply(data, 1, makeplot)

#cat("Finish plotting... ", file=logFile, append=TRUE, sep = "\n")

#testdata2<-td[[1]]
#print("Done!!!")

#export sys
#Recomended to save the session info for reproducibility details
#Create a file sessionInfo.txt and save the R sessionInfo()
writeLines(capture.output(sessionInfo()), paste(outDir,.Platform$file.sep,"sessionInfo.txt",sep = ""))

