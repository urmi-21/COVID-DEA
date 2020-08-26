library(readr)
library(dplyr)
library(plyr)
library(reshape)
library(data.table)



metadataFile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_MetaData_7142_23_updated.tsv"
infile<-"C:/Users/mrbai/Desktop/mog_demo/tr2/mogData_R_Data.txt"
outDir<-"C:/Users/mrbai/Desktop/mog_demo/anova_out_full/"

#read data
data<-read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
allGenes<-data$Name
#read metadata
metadata<- read_delim(metadataFile, "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(metadata)[which(colnames(metadata)=="portions.analytes.aliquots.submitter_id")]='SID'

#use only
races<-c("black or african american","white")
#filter metadata to keep only data columns and races
metadata<-metadata%>%filter(SID %in% colnames(data)) %>% filter(clinical.race %in% races)
#change metadata disease column
metadata<-metadata%>%mutate(clinical.disease=gsub("_normal", "", clinical.disease))
diseases<-as.character(unique(metadata$clinical.disease))

#metadata categories
metadata <- metadata %>% mutate(age_new=ifelse(age<=25,'A',ifelse(age<=50,'B',ifelse(age<=75,'C','D'))))

#t<-metadata%>%select(age,age_new)

removeNARows <- function(data, colsToCheck) {
  result <- complete.cases(data[, colsToCheck])
  return(data[result, ])
}


anova_rowwise <- function(data_row,testCols){
  #do log transform
  data_row<-log(data_row+1)
  genemelt<-melt(data_row)
  genemelt$SID<-rownames(genemelt)
  
  #add metadata
  genemelt <- left_join(genemelt,metadata,by="SID")
  
  #lapply(genemelt, levels)
  #return(genemelt)
  #do anova
  #testCols<-c("clinical.race","clinical.gender","clinical.disease","age_new")
  #remove any rows with NA values
  genemelt<-removeNARows(genemelt,testCols)
  #formula
  f <- as.formula(paste('value',paste(testCols, collapse = " + "),sep = " ~ "))
  #print(f)
  fit <- lm(f, data=genemelt)
  anova(fit)
  #fit$coefficients
}

anova_listwise <- function(df,tcols) {
  #get data for samples in df; 1st col is gene name
  print(unique(df$clinical.disease))
  thisData<-data[,c('Name',df$SID)]
  genesName<-thisData$Name
  thisData$Name<-NULL
  out<-apply(thisData, 1, anova_rowwise,testCols=tcols)
  names(out)<-genesName
  return(out)
}

write_list_results <- function(fname){
  print(paste('in',fname))
  thisRes<-lapp_out[[fname]]
  finalres<-data.frame(Gene=character(),
                       pv_race=double(), 
                       pv_gender=double(),
                       pv_age=double(),
                       Residual=double(),
                       stringsAsFactors=FALSE) 
  
  for(i in 1:length(thisRes)){
    generes=thisRes[[i]]
    thisGene<-names(thisRes)[i]
    #print(generes)
    pv=generes$`Pr(>F)`
    ss=generes$`Sum Sq`
    finalres[nrow(finalres) + 1,] = c(thisGene,pv[1],pv[2],pv[3],ss[4])
  }
  #write to file
  print(paste('writing',fname))
  write_tsv(finalres,paste(outDir,.Platform$file.sep,fname,"_AnovaOut.tsv",sep = ""))
}

write_list_results2 <- function(fname){
  print(paste('in',fname))
  thisRes<-lapp_out[[fname]]
  finalres<-data.frame(Gene=character(),
                       pv_race=double(), 
                       pv_gender=double(),
                       pv_condition=double(),
                       pv_age=double(),
                       Residual=double(),
                       stringsAsFactors=FALSE) 
  
  for(i in 1:length(thisRes)){
    generes=thisRes[[i]]
    thisGene<-names(thisRes)[i]
    #print(generes)
    pv=generes$`Pr(>F)`
    ss=generes$`Sum Sq`
    finalres[nrow(finalres) + 1,] = c(thisGene,pv[1],pv[2],pv[3],pv[4],ss[5])
  }
  #write to file
  print(paste('writing',fname))
  write_tsv(finalres,paste(outDir,.Platform$file.sep,fname,"_AnovaOut.tsv",sep = ""))
}

#########################################################################################

#split data by condition
#exclude breast,uterus,prostate
metadata<-metadata%>%filter(! (clinical.disease %in% c("Breast","Uterus","Prostate","BRCA","UCEC","UCS","PRAD")))
splitset<-split(metadata,list(metadata$clinical.disease))
#test<-list(splitset[[1]],splitset[[1]])
testCols<-c("clinical.race","clinical.gender","age_new")
lapp_out<-lapply(splitset, anova_listwise,tcols=testCols)
#write to file
sapply(names(lapp_out), write_list_results)


#do for gtex,tcga
#include only selected disease
tokeep<-c("BRCA","UCEC","COAD","THCA","LUAD","LUSC","KIRC","KIRP","Breast","Thyroid","Lung","Esophagus","Stomach","Colon","Prostate","Liver","Kidney" ,"Uterus")
metadata<-metadata%>%filter(clinical.disease %in% tokeep)
#split data by condition
splitset<-split(metadata,list(metadata$Project))
#test<-list(splitset[[1]],splitset[[1]])
testCols<-c("clinical.race","clinical.gender","clinical.disease","age_new")
lapp_out<-lapply(splitset, anova_listwise,tcols=testCols)
#write to file
sapply(names(lapp_out), write_list_results2)

#########################################################################################
