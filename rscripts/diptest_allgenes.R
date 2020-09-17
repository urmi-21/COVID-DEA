library(readr)
library(dplyr)
library(diptest)
library(plyr)
library(data.table)



#This script takes 3 arguments
#arg1: Path to input file (this file is created by MOG)
#arg2: Path to the metadatafile (used in the MOG project)
#arg3: Path to the output directory (specified by the MOG user in MOG)

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Invalid Arguments", call.=FALSE)
} 

#read data
#infile<-args[1]
#metadataFile<-args[2]
#outDir<-args[3]



metadataFile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_MetaData_7142_23_updated.tsv"
infile<-"C:/Users/mrbai/Desktop/mog_demo/tr2/mogData_R_Data.txt"
outDir<-"C:/Users/mrbai/Desktop/mog_demo/tr2/"

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
metadata<-metadata %>% filter(clinical.race %in% races)

#change metadata disease column
metadata<-metadata%>%mutate(clinical.disease=gsub("_normal", "", clinical.disease))
diseases<-as.character(unique(metadata$clinical.disease))


#############Start###############
dodiptest <- function(df) {
  #get data for samples in df; ist col is gene name
  thisData<-data[,c('Name',df$SID)]
  dt_out<-apply(thisData, 1, diptest_rowwise)
  dt_out<-t(data.frame(dt_out))
  
}

diptest_rowwise <- function(data_row){
  #print(data_row[1])
  thisname<-data_row[1]
  thisdata<-data_row[2:length(data_row)]
  tmp<-dip.test(log2(as.numeric(thisdata)+1))
  c(Gene=thisname,N=tmp$nobs,Method=tmp$method,Statistic=tmp$statistic,pvalue=format(round(tmp$p.value, 4), nsmall = 4))
}

#splitset<-split(metadata,list(metadata$clinical.race,metadata$clinical.disease))
splitset<-c(split(metadata,list(metadata$clinical.race,metadata$Project)),split(metadata,list(metadata$clinical.race,metadata$clinical.disease)))
lapp_out<-lapply(splitset, dodiptest)


#write dip test res to file
for (n in 1:length(lapp_out)){
  #print(n)
  thisname<-names(lapp_out)[n]
  xd<-data.frame(lapp_out[[n]],stringsAsFactors = F)
  xd<-as.data.frame(xd)
  temp<-strsplit(thisname, '[.]')[[1]]
  xd$Race<-temp[1]
  xd$Condition<-temp[2]
  xd$padj<-p.adjust(xd$pvalue,method="fdr")
  #write to file
  write_tsv(xd,paste(outDir,temp[1],"_",temp[2],"_Diptest.tsv",sep=""))
}



###KS test####
dokstest <- function(df) {
  #get data for samples in df; ist col is gene name
  thisData<-data[,c('Name',df$SID)] #all genes in this disease
  blackSamps<-df %>% filter(clinical.race == "black or african american") %>% select(SID)
  whiteSamps<-df %>% filter(clinical.race == "white") %>% select(SID)
  
  dt_out<-apply(thisData, 1, kstest_rowwise,b=blackSamps,w=whiteSamps)
  dt_out<-t(data.frame(dt_out))
  #newdf<-data.frame(c(Race=unique(df$clinical.race),Disease=unique(df$clinical.disease)))
  #print(newdf)
  #data.frame(dt_out,newdf)
  
}

kstest_rowwise <- function(data_row,b,w){
  #print(data_row[1])
  thisname<-data_row[1]
  thisdata<-data_row[2:length(data_row)]
  blackData<-thisdata[b$SID]
  whiteData<-thisdata[w$SID]
  #ks.test
  tmp<-ks.test(log2(as.numeric(blackData)+1),log2(as.numeric(whiteData)+1))
  c(Gene=thisname,N_black=length(blackData),N_white=length(whiteData),Method=tmp$method,Statistic=tmp$statistic,pvalue=format(round(tmp$p.value, 4), nsmall = 4))
}

#splitset<-split(metadata,list(metadata$clinical.disease))
splitset<-c(split(metadata,list(metadata$Project)),split(metadata,list(metadata$clinical.disease)))
lapp_out<-lapply(splitset, dokstest)
#write ks test res to file
for (n in 1:length(lapp_out)){
  #print(n)
  thisname<-names(lapp_out)[n]
  xd<-data.frame(lapp_out[[n]],stringsAsFactors = F)
  xd<-as.data.frame(xd)
  xd$Condition<-thisname
  xd$padj<-p.adjust(xd$pvalue,method="fdr")
  #write to file
  write_tsv(xd,paste(outDir,thisname,"_KStest",".tsv",sep=""))
}


#Recomended to save the session info for reproducibility details
#Create a file sessionInfo.txt and save the R sessionInfo()
writeLines(capture.output(sessionInfo()), paste(outDir,.Platform$file.sep,"sessionInfo.txt",sep = ""))


