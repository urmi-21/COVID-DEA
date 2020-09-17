library(readr)
library(dplyr)
library(diptest)
library(plyr)
library(scales)
library(data.table)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(grid)  
library(gtable)  



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



#metadataFile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_MetaData_7142_23_updated.tsv"
#infile<-"C:/Users/mrbai/Desktop/mog_demo/tcga_gt_vp/mogData_R_Data.txt"
#outDir<-"C:/Users/mrbai/Desktop/mog_demo/tcga_gt_vp/"

#read data
data<-read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
allGenes<-data$Name
#read metadata
metadata<- read_delim(metadataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(metadata)[which(colnames(metadata)=="portions.analytes.aliquots.submitter_id")]='SID'

#filter metadata to keep only data columns
metadata<-metadata%>%filter(SID %in% colnames(data))

#races<-as.character(unique(metadata$clinical.race))
#use only
races<-c("black or african american","white")

#change metadata disease column
metadata<-metadata%>%mutate(clinical.disease=gsub("_normal", "", clinical.disease))
diseases<-as.character(unique(metadata$clinical.disease))


#create all gtex data
gtexSamps<-metadata %>% filter(Project == 'GTEX') %>% filter(SID %in% colnames(data))
gtexData<-data%>%select("Name",gtexSamps$SID)


#create all tcga data
tcgaSamps<-metadata %>% filter(Project == 'TCGA') %>% filter(SID %in% colnames(data))
tcgaData<-data%>%select("Name",tcgaSamps$SID)


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
#Gtex
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','CCL3L3',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','CCL3L3',1.51E-04)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','CCL3L3',0.004957488)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','CCL3L3',0.902384316)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','CCL3L3',0.267891803)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','CCL3L3',0.010467593)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','CCL3L3',0.219733335)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','CCL3L3',0.519103523)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','CCL3L3',0.014593547)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','CCL3L3',0.815962134)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','CXCL9',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','CXCL9',0.842110806)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','CXCL9',0.944375761)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','CXCL9',0.822358748)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','CXCL9',0.810958813)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','CXCL9',0.950262726)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','CXCL9',0.293670122)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','CXCL9',0.993993027)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','CXCL9',0.879217924)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','CXCL9',0.872817092)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','CXCL10',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','CXCL10',0.854369476)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','CXCL10',0.959059035)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','CXCL10',0.832170336)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','CXCL10',0.808801749)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','CXCL10',0.695339572)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','CXCL10',0.481972356)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','CXCL10',0.985476319)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','CXCL10',0.699985299)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','CXCL10',0.932729814)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','F8A1',0.214107552)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','F8A1',1.43E-05)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','F8A1',6.72E-11)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','F8A1',0.63686807)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','F8A1',0.165602522)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','F8A1',0.001150247)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','F8A1',0.202835382)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','F8A1',0.002191171)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','F8A1',8.67E-07)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','F8A1',0.68532559)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','F8A2',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','F8A2',0.0029254)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','F8A2',2.02E-05)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','F8A2',0.634945875)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','F8A2',0.069790598)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','F8A2',7.83E-04)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','F8A2',0.219733335)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','F8A2',0.215134408)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','F8A2',3.01E-04)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','F8A2',0.757166559)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','F8A3',0.414838889)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','F8A3',8.74E-04)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','F8A3',2.07E-10)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','F8A3',0.737811502)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','F8A3',0.509658504)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','F8A3',2.00E-05)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','F8A3',0.352800368)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','F8A3',0.012826583)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','F8A3',1.57E-06)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','F8A3',0.670802176)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','GSTM1',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','GSTM1',0.681084103)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','GSTM1',0.014995615)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','GSTM1',0.912382704)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','GSTM1',0.484218347)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','GSTM1',0.001148033)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','GSTM1',0.680068434)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','GSTM1',0.874012114)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','GSTM1',0.011825393)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','GSTM1',0.53232376)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','IL6ST',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','IL6ST',0.904753559)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','IL6ST',0.948109124)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','IL6ST',0.923971046)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','IL6ST',0.976761683)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','IL6ST',0.15781366)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','IL6ST',0.970755751)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','IL6ST',0.985476319)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','IL6ST',0.988133044)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','IL6ST',0.92407056)
MOG_MW[nrow(MOG_MW) + 1,] =c('Breast','TMPRSS6',0.996800546)
MOG_MW[nrow(MOG_MW) + 1,] =c('Colon','TMPRSS6',0.995665763)
MOG_MW[nrow(MOG_MW) + 1,] =c('Esophagus','TMPRSS6',0.999652187)
MOG_MW[nrow(MOG_MW) + 1,] =c('Kidney','TMPRSS6',0.9544692)
MOG_MW[nrow(MOG_MW) + 1,] =c('Liver','TMPRSS6',0.995648507)
MOG_MW[nrow(MOG_MW) + 1,] =c('Lung','TMPRSS6',0.640421868)
MOG_MW[nrow(MOG_MW) + 1,] =c('Prostate','TMPRSS6',0.464519136)
MOG_MW[nrow(MOG_MW) + 1,] =c('Stomach','TMPRSS6',0.992248175)
MOG_MW[nrow(MOG_MW) + 1,] =c('Thyroid','TMPRSS6',0.64091724)
MOG_MW[nrow(MOG_MW) + 1,] =c('Uterus','TMPRSS6',0.752606253)

#All
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','CCL3L3',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','CCL3L3',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','CCL3L3',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','CXCL9',0.075808844)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','CXCL9',0.525008244)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','CXCL9',0.001334638)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','CXCL10',0.696101602)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','CXCL10',0.541004307)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','CXCL10',0.132690647)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','F8A1',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','F8A1',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','F8A1',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','F8A2',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','F8A2',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','F8A2',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','F8A3',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','F8A3',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','F8A3',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','GSTM1',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','GSTM1',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','GSTM1',1.60E-13)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','IL6ST',2.29E-12)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','IL6ST',0.56214713)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','IL6ST',0)
MOG_MW[nrow(MOG_MW) + 1,] =c('ALL','TMPRSS6',0.029311269)
MOG_MW[nrow(MOG_MW) + 1,] =c('GTEx','TMPRSS6',0.877484586)
MOG_MW[nrow(MOG_MW) + 1,] =c('TCGA','TMPRSS6',0.008960054)





MOG_MW$FDR<-as.numeric(MOG_MW$FDR)

#cat("Finish loading MW results... ", file=logFile, append=TRUE, sep = "\n")

#create plots
#do for each row/gene in data
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
  
  
  #add gtex all by copying all data
  tissuelist<-c("Breast","Colon","Esophagus","Kidney","Liver","Lung","Prostate","Stomach","Thyroid","Uterus")
  tumorlist<-c("BRCA","COAD","KIRC","KIRP","LUAD","LUSC","THCA","UCEC")
  
  testdata2gtex<-testdata2 %>% filter(Disease %in% tissuelist)
  testdata2tcga<-testdata2 %>% filter(Disease %in% tumorlist)
  testdata2gtex<-testdata2gtex%>%mutate(Disease="GTEx")
  testdata2tcga<-testdata2tcga%>%mutate(Disease="TCGA")
  
  
  #infer gtex or tcga
  #allDis<-as.character(unique(testdata2$Disease))
  #quick test
  #if('Colon' %in% allDis){
  #  allName='GTEx'
  #  xtitle="Tissue"
  #}else{
  #  allName='TCGA'
  #  xtitle="Tumor"
  #}
  #testdata2copy<-testdata2copy%>%mutate(Disease=allName)
  
  #calculate FC for 'pooled' samples
  ##GTEx
  thisGene=genename
  grpaa<-testdata2gtex%>%filter(Race=='AA')
  grpea<-testdata2gtex%>%filter(Race=='EA')
  aadata=as.matrix(grpaa[2])
  eadata=as.matrix(grpea[2])
  maa<-mean(aadata)
  mea<-mean(eadata)
  logfc=mean(log2(aadata+1))-mean(log2(eadata+1))
  fc=2^logfc
  FCResults[nrow(FCResults) + 1,] = c('GTEx',thisGene,length(grpaa[2]),length(grpea[2]),maa,mea,logfc,fc)
  #perform diptest aa
  tmp<-dip.test(log2(aadata+1))
  dipTestResults[nrow(dipTestResults) + 1,] = c('GTEx','black or african american',thisGene,tmp$data.name,tmp$nobs,tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
  #perform diptest ea
  tmp<-dip.test(log2(eadata+1))
  dipTestResults[nrow(dipTestResults) + 1,] = c('GTEx','white',thisGene,tmp$data.name,tmp$nobs,tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
  #perform KS
  tmp<-ks.test(log2(aadata+1),log2(eadata+1))
  dim(KSTestResults)
  KSTestResults[nrow(KSTestResults) + 1,] = c('GTEx',thisGene,tmp$data.name,length(blackData),length(whiteData),tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
  dim(KSTestResults)
  ##TCGA
  thisGene=genename
  grpaa<-testdata2tcga%>%filter(Race=='AA')
  grpea<-testdata2tcga%>%filter(Race=='EA')
  aadata=as.matrix(grpaa[2])
  eadata=as.matrix(grpea[2])
  maa<-mean(aadata)
  mea<-mean(eadata)
  logfc=mean(log2(aadata+1))-mean(log2(eadata+1))
  fc=2^logfc
  FCResults[nrow(FCResults) + 1,] = c('TCGA',thisGene,length(grpaa[2]),length(grpea[2]),maa,mea,logfc,fc)
  #perform diptest aa
  tmp<-dip.test(log2(aadata+1))
  dipTestResults[nrow(dipTestResults) + 1,] = c('TCGA','black or african american',thisGene,tmp$data.name,tmp$nobs,tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
  #perform diptest ea
  tmp<-dip.test(log2(eadata+1))
  dipTestResults[nrow(dipTestResults) + 1,] = c('TCGA','white',thisGene,tmp$data.name,tmp$nobs,tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
  #perform KS
  tmp<-ks.test(log2(aadata+1),log2(eadata+1))
  dim(KSTestResults)
  KSTestResults[nrow(KSTestResults) + 1,] = c('TCGA',thisGene,tmp$data.name,length(blackData),length(whiteData),tmp$method,tmp$statistic,format(round(tmp$p.value, 4), nsmall = 4))
  dim(KSTestResults)
  
  #update global vars
  KSTestResults<<-KSTestResults
  FCResults<<-FCResults
  dipTestResults<<-dipTestResults
  
  #combine pooled data in diseasewise data
  testdata2<-bind_rows(testdata2,testdata2tcga)
  testdata2<-bind_rows(testdata2,testdata2gtex)
  #return(testdata2)
  
  
  
  medians<-aggregate( get(genename) ~ Disease+Race, testdata2, median )
  colnames(medians)[3]<-'median'
  means<-aggregate( get(genename) ~ Disease+Race, testdata2, mean )
  colnames(means)[3]<-'mean'
  
  #plot violin plot
  dodge <- position_dodge(width = 0.9)
  
  ###Do calculations to add significance
  testdata2$Race<- as.factor(testdata2$Race)
  testdata2$Disease<-as.factor(testdata2$Disease)
  
  tdst<-compare_means(
    as.formula(paste(genename,"~ Race")), data = testdata2, group.by = "Disease",
    method = "t.test", ref.group = "AA"
  )
  
  
  #max val used for position to put significance results
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
  FCval$FC<-format(round(as.numeric(FCval$FC), 1), nsmall = 1)# round the value to one decimal
  #FCval$FC<-ceiling(as.numeric(FCval$FC)) #to integer
  #FCval$FC<-round_any(as.numeric(FCval$FC), 0.1) #to nearest 0.1
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
  
  
  #custom order
  #dlevels<-levels(testdata2$Disease)
  #print(dlevels)
  newLevels<-c("GTEx","Breast","Colon","Esophagus","Kidney","Liver","Lung","Prostate","Stomach","Thyroid","Uterus",
               "TCGA","BRCA","COAD","KIRC","KIRP","LUAD","LUSC","THCA","UCEC"
  )
  
  #newLevels[length(newLevels)+1]<-allName
  #print(newLevels)
  testdata2$Disease <- factor(testdata2$Disease, levels = newLevels)
  #return(testdatestdata2$Diseaseta2)
  means$Disease<- factor(means$Disease, levels = newLevels)
  #print(means)
  #return(testdata2)
  testdata2[,genename]<-testdata2[,genename]+1
  plot1<-ggplot(testdata2, aes_string(x='Disease',  y=genename))+
    geom_violin(aes(fill=Race),position = dodge,alpha = 0.1,show.legend = FALSE)+
    geom_boxplot(aes(fill=Race),width=.1, outlier.colour=NA, position = dodge,show.legend = FALSE)+
    stat_pvalue_manual(dip,x='Disease',color = "Blue",label = "p.signif",size = 7)+
    #stat_pvalue_manual(KS,x='Disease',color = "firebrick3",label = "p.signif",size = 7)+
    stat_pvalue_manual(MW,x='Disease',color = "seagreen",label = "p.signif",size = 7)+
    stat_pvalue_manual(FCval,x='Disease',color = "Black",label = "FC",size = 4.5)+
    scale_y_log10()+
    theme(legend.key=element_blank(),legend.position= c(0.92, 0.35),
          #legend.spacing.x = unit(0.5, 'cm'),#legend space
          legend.title=element_blank(),legend.text = element_text(size=10,margin = margin(t=1,r=40,b=1,l=1)), #top legend
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(size = 10,face = "bold",angle = 0, vjust = 0.5, hjust=1), #x axis labels
          axis.text.y = element_text(size = 10,face = "bold",angle = 0, vjust = 0.5, hjust=1), #y axis labels
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ #the y axis label
    ylab(paste(genename,"(FPKM)"))+xlab("")+
    geom_hline(data = means, aes(yintercept = mean, color = Race),size=1,show.legend = FALSE)+ 
    guides(colour=FALSE)+guides(color = guide_legend(override.aes = list(size = 3)))+
    #facet_wrap( ~ Disease, scales="free_x",ncol = 12)+ for same y axis
    facet_wrap( ~ Disease, scales="free",ncol = 11)+
    xlab("")+
    theme(strip.background = element_blank(),strip.text.x = element_blank(),
          panel.spacing = unit(0.1, "lines"),panel.margin.x=unit(0.1, "lines"),
          panel.spacing.y=unit(1, "lines")
    )
  #return(plot1)
  
  #customlegend
  # Construct the six grobs - three symbols and three labels
  L0 = rectGrob(height = .5, width = .5, gp = gpar(fill = "#F8766D", col = NA))
  L1 = rectGrob(height = .5, width = .5, gp = gpar(fill = "#00BFC4", col = NA))
  L4 = textGrob("*", x = .2,y=0.3, just = "left",gp=gpar(fontsize = 22,col="blue",fontface = 'bold'))
  L2 = textGrob("*", x = .2,y=0.3, just = "left",gp=gpar(fontsize = 22,col="darkgreen",fontface = 'bold'))
  
  T0 = textGrob("AA", x = 0.05, just = "left",gp=gpar(fontsize = 11,col="black",fontface = 'bold'))
  T1 = textGrob("EA", x = .05, just = "left",gp=gpar(fontsize = 11,col="black",fontface = 'bold'))
  T4 = textGrob("Dip test significant", x = .05, just = "left",gp=gpar(fontsize = 11,col="blue",fontface = 'bold'))
  T2 = textGrob("MW test significant", x = .05, just = "left",gp=gpar(fontsize = 11,col="darkgreen",fontface = 'bold'))
  
  # Construct a gtable - 2 columns X 4 rows
  leg = gtable(width = unit(c(0.5,4), "cm"), height = unit(c(1,0.5,0.5,0.5,0.5,0.5), "cm"))
  leg = gtable_add_grob(leg, rectGrob(gp = gpar(fill = NA, col = "black")), t=2,l=1,b=6,r=2)
  
  # Place the six grob into the table
  leg = gtable_add_grob(leg, L0, t=2, l=1)
  leg = gtable_add_grob(leg, L1, t=3, l=1)
  leg = gtable_add_grob(leg, L2, t=5, l=1)
  leg = gtable_add_grob(leg, L4, t=6, l=1)
  
  leg = gtable_add_grob(leg, T0, t=2, l=2)
  leg = gtable_add_grob(leg, T1, t=3, l=2)
  leg = gtable_add_grob(leg, T2, t=5, l=2)
  leg = gtable_add_grob(leg, T4, t=6, l=2)
  
  # Give it a title (if needed)
  leg = gtable_add_grob(leg, textGrob("Legend"), t=1, l=1, r=2)
  # Get the ggplot grob for plot1
  g = ggplotGrob(plot1)
  # Get the position of the panel,
  # add a column to the right of the panel, 
  # put the legend into that column, 
  # and then add another spacing column
  pos = g$layout[grepl("panel-11-2", g$layout$name), c('t', 'l')]
  pos$t<-13.2
  pos$l<-43.8
  pos$z<-1
  #g = gtable_add_cols(g, sum(leg$widths), pos$l)
  g = gtable_add_grob(g, leg, t = pos$t, l = pos$l, z = Inf )
  # Draw it
  #grid.newpage()
  #grid.draw(g)
  
  #print(paste('saving:',genename))
  #plot
  ggsave(paste(outDir,.Platform$file.sep,genename,".png",sep = ""),g,width = 13, height = 7, units = "in",dpi = 300)
  
}

apply(data, 1, makeplot)

#save results
write_tsv(dipTestResults,paste(outDir,.Platform$file.sep,'diptest.tsv',sep = ""))
write_tsv(KSTestResults,paste(outDir,.Platform$file.sep,'KStest.tsv',sep = ""))
write_tsv(FCResults,paste(outDir,.Platform$file.sep,'FCResults.tsv',sep = ""))

#Recomended to save the session info for reproducibility details
#Create a file sessionInfo.txt and save the R sessionInfo()
writeLines(capture.output(sessionInfo()), paste(outDir,.Platform$file.sep,"sessionInfo.txt",sep = ""))


