library(readr)
library(plyr)
library(dplyr)
library(reshape)
library(data.table)
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(tidyr)




#read data; download data from http://metnetweb.gdcb.iastate.edu/MetNet_MetaOmGraph_download.php
metadataFile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_MetaData_7142_23_updated.tsv"
infile<-"D:/MOGdata/mog_testdata/cancer/raw/normalized/alldata_reduced/TCGA_GTEX_MOG/TCGA_GTEX_MOG_proj/US_5-20-20_GTEx_TCGA_MOG/TCGA_GTEX_Data_18212_7142.tsv"
outDir<-"C:/Users/mrbai/Desktop/mog_demo/limma_BRCAsubtypes/"

data<-read_delim(infile,"\t", escape_double = FALSE, trim_ws = TRUE)
colnames(data)[1]<-"Name"
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


#######Get subtypes#########
subtypes <- PanCancerAtlas_subtypes()

subtypes_tab<-DT::datatable(subtypes,
                            filter = 'top',
                            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
                            rownames = FALSE)
subtypes_tab<-as.data.frame(subtypes)

colnames(subtypes_tab)[1]<-"SID"
metadata<-left_join(metadata,subtypes_tab,by="SID")
################################
#keep only brca with subtypes
metadata <- metadata %>% filter(clinical.disease %in% c('BRCA')) %>%
  filter(Subtype_Selected %in% c("BRCA.LumA","BRCA.Her2","BRCA.LumB","BRCA.Basal"))
#remove NA rows in factor of interest
metadata<-metadata %>% drop_na(clinical.race, age_new,Subtype_Selected)



#keep all samples in metadata and in order
data<-data[,c("Name",metadata$SID)]
genes<-data[,1]
data<-data[,2:ncol(data)]

#check order of samples and factors
if (!all(colnames(data)==metadata$SID)){
  print("ERROR in data order")
}


logData<-log(data+1)

#create design
race<-as.factor(metadata$clinical.race)
age<-factor(metadata$age_new)
subtype<-factor(metadata$Subtype_Selected)


#design<-model.matrix(~race+age+subtype)
design<-model.matrix(~0+race+age+subtype)
colnames(design)[1:2]<-c("AA","EA")
cont.mat <- makeContrasts(AA-EA, levels = design)

#fit lm
fit <- lmFit(logData, design)
#efit <- eBayes(fit, trend=TRUE,robust=TRUE)
#tt<-topTable(efit, coef=ncol(design),number=50)
#tt$gene<-genes[as.numeric(rownames(tt)),1]$Name

cont.fit <- contrasts.fit(fit, cont.mat)
cont.fit.eb <- eBayes(cont.fit, trend = TRUE)
cont.res <- topTable(cont.fit.eb,coef=1,n = Inf, adjust="fdr")
cont.res$Name<-genes[as.numeric(rownames(cont.res)),1]$Name
#sort by FC
cont.res <- cont.res %>% arrange(desc(logFC))
#reorder columns
cont.res<-cont.res%>%select(c(Name,AveExpr,t,logFC, P.Value,    adj.P.Val, B))
#rename columns
colnames(cont.res)[6]<-"Adj pval"


#save to file
write_tsv(cont.res,paste(outDir,.Platform$file.sep,"BRCA_limmaOut.tsv",sep = ""))
