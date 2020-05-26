#!/usr/bin/env Rscript

#author: urmi
#A script to plot violin plots. Invoked through MetaOmGraph.

library(dplyr)
library(plyr)
library(scales)
library(readr)
library(data.table)
library(ggplot2)
library(ggthemes)

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

#infile<-'covidset1.txt'
#metadataFile<-'D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/TCGA_GTEX_MetaData_7142_23.tsv'
#outDir<-'D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/hmout11'

#read metadata
metadata<- read_delim(metadataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(metadata)[which(colnames(metadata)=="portions.analytes.aliquots.submitter_id")]='SID'

#read data
data<-read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
#rownames(data)<-data$Name
#data$Name<-NULL
#colnames(data)[which(colnames(data)=="Name")]='SID'


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
  print(paste('dim:',dim(tdata),' gene:',genename))
  
  #join with metadata
  tdata<-merge(tdata,metadata,all.x = TRUE,by='SID')
  #return(tdata)
  testdata2<- tdata %>% select('clinical.disease',genename,'clinical.race')
  colnames(testdata2)[1]='Disease'
  colnames(testdata2)[3]='Race'
  
  print(paste('dim:',dim(testdata2)))
  
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
  
  
  
  
  #plot violin plot
  dodge <- position_dodge(width = 0.9)
  
  ggplot(testdata2, aes_string(x='Disease',  y=genename ,fill='Race'))  +
    geom_violin(position = dodge,alpha = 0.1,show.legend = FALSE)+
    geom_boxplot(width=.1, outlier.colour=NA, position = dodge,fatten = NULL,show.legend = FALSE)+
    #geom_jitter(position=position_jitter(0.2),aes_string(color = 'Race'))+
    scale_y_log10()+
    theme(legend.key=element_blank(),legend.position="top",legend.title=element_blank(),legend.text = element_text(size=7),axis.ticks.x=element_blank(),axis.text.x = element_text(size = 6,face = "bold"),axis.text.y = element_text(size = 6,face = "bold"),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), axis.title=element_text(size=7,face="bold"))+
    
    ylab(paste(genename,"(FPKM)"))+xlab("")+
    #coord_cartesian(ylim = c(0, 500))+
    geom_hline(data = means, aes(yintercept = mean, color = Race,show.legend = FALSE))+ 
    guides(colour=FALSE)+guides(color = guide_legend(override.aes = list(size = 3)))+
    #geom_hline(data = medians, aes(yintercept = median))+
    facet_wrap( ~ Disease, scales="free",ncol = 8)+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),panel.spacing = unit(0.1, "lines"),panel.margin.x=unit(0.1, "lines")
    )
  
  
  print(paste('saving:',genename))
  #plot
  ggsave(paste(outDir,.Platform$file.sep,genename,".png",sep = ""),width = 8, height = 3, units = "in",dpi = 300)
  
}

apply(data, 1, makeplot)
print("Done!!!")





