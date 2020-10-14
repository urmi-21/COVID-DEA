#Perform DE analysis using limma tissue/tumor wise

library(readr)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape)
library(data.table)
library(limma)
library(edgeR)




metadataFile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_MetaData_7142_23_updated_bmi.tsv"
infile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_Data_18212_7142.tsv"
outDir<-"C:/Users/mrbai/Desktop/mog_demo/limma_tissuewise_bmi2/"
#create outdir
dir.create(file.path(outDir))

#read data
data<-read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
colnames(data)[1]<-"Name"
allGenes<-data$Name
#read metadata
metadata<- read_delim(metadataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
#change sample ID column name to SID
colnames(metadata)[which(colnames(metadata)=="portions.analytes.aliquots.submitter_id")]='SID'
#Filter samples by race
races<-c("black or african american","white")
#filter metadata to filter samples by Race
metadata<-metadata %>% filter(clinical.race %in% races)
#rename race categories
metadata<-metadata %>% mutate(clinical.race=ifelse(clinical.race=="white","EA","AA"))
#change metadata disease column
metadata<-metadata%>%mutate(clinical.disease=gsub("_normal", "", clinical.disease))
diseases<-as.character(unique(metadata$clinical.disease))
#create age factor
metadata <- metadata %>% mutate(age_new=ifelse(age<=25,'A',ifelse(age<=50,'B',ifelse(age<=75,'C','D'))))
#remove NA rows from factors of interest
metadata<-metadata %>% drop_na(clinical.race,age_new,clinical.gender,BMIcat)

dplyr::count(metadata,clinical.race,clinical.disease)


##############################Funcs###############################################

limma_listwise <- function(df,cofacs) {
  print(unique(df$clinical.disease))
  
  #get data; cols are in order of metadata
  thisData<-data[,c('Name',df$SID)]
  genes<-thisData[,1]
  thisData<-thisData[,2:ncol(thisData)]
  #check order matches with metadata
  if (!all(colnames(thisData)==df$SID)){
    print("ERROR in data order")
    return(NULL)
  }
  
  #log data
  logData<-log(thisData+1)
  
  #create design
  #race is primary factor
  race<-as.factor(df$clinical.race)
  
  #other factors; create vectors for cofacs
  for(x in cofacs){
    #print(x)
    assign(x,as.factor(df[[x]]))
  }
  
  f <- as.formula(paste0('~0+race+',paste(cofacs, collapse = "+")))
  
  design<-model.matrix(f)
  colnames(design)[1:2]<-c("AA","EA")
  cont.mat <- makeContrasts(AA-EA, levels = design)
  
  #start limma
  #fit lm
  fit <- lmFit(logData, design)
  cont.fit <- contrasts.fit(fit, cont.mat)
  cont.fit.eb <- eBayes(cont.fit, trend = TRUE, robust=TRUE)
  cont.res <- topTable(cont.fit.eb,coef=1,n = Inf, adjust="fdr")
  cont.res$Name<-genes[as.numeric(rownames(cont.res)),1]$Name
  #sort by FC
  cont.res <- cont.res %>% arrange(desc(logFC))
  #reorder columns
  cont.res<-cont.res%>%select(c(Name,AveExpr,t,logFC, P.Value,    adj.P.Val, B))
  #rename columns
  colnames(cont.res)[6]<-"Adj pval"
  #save to file
  #write_tsv(cont.res,paste(outDir,.Platform$file.sep,"BRCA_limmaOut.tsv",sep = ""))
  
  return(cont.res)
}


write_list_results <- function(fname){
  thisRes<-lapp_out[[fname]]
  print(paste('writing',fname))
  write_tsv(thisRes,paste(outDir,.Platform$file.sep,fname,"_limmaOut.tsv",sep = ""))
}


#########################################################################################

#split data by condition
#exclude breast,uterus,prostate
metadata_ss<-metadata%>%filter(! (clinical.disease %in% c("Breast","Uterus","Prostate","BRCA","UCEC","UCS","PRAD")))
splitset<-split(metadata_ss,list(metadata_ss$clinical.disease))
#specify factors other than race
testCols<-c("clinical.gender","age_new","BMIcat")
lapp_out<-lapply(splitset, limma_listwise,cofacs=testCols)
#write to file
sapply(names(lapp_out), write_list_results)

#do for gender specific samples using race and age
metadata_gs<-metadata%>%filter((clinical.disease %in% c("Breast","Uterus","Prostate","BRCA","UCEC","UCS","PRAD")))
splitset<-split(metadata_gs,list(metadata_gs$clinical.disease))
#specify factors other than race
testCols<-c("age_new","BMIcat")
lapp_out<-lapply(splitset, limma_listwise,cofacs=testCols)
#write to file
sapply(names(lapp_out), write_list_results)




#do for gtex,tcga
#include only selected disease
tokeep<-c("BRCA","UCEC","COAD","THCA","LUAD","LUSC","KIRC","KIRP","Breast","Thyroid","Lung","Esophagus","Stomach","Colon","Prostate","Liver","Kidney" ,"Uterus")
metadata_tg<-metadata%>%filter(clinical.disease %in% tokeep)
#split data by condition
splitset<-split(metadata_tg,list(metadata_tg$Project))
#test<-list(splitset[[1]],splitset[[1]])
testCols<-c("clinical.gender","clinical.disease","age_new","BMIcat")
lapp_out<-lapply(splitset, limma_listwise,cofacs=testCols)
#write to file
sapply(names(lapp_out), write_list_results)

#########################################################################################

