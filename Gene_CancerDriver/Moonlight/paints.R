###Moonlight Cancer  driver genes lung cancer subtype luad, subsubtype luad1

#Import Moonlight
library(MoonlightR)
library(tidyverse)
library(ggplot2)

#data load
dataDEGs <- get(load("../inputfiles/LUAD_dataDEGs_deconvolution.rda")) 

DataFilt <- get(load("../inputfiles/LUAD_dataFilt_adjusted.rda"))

#Find normal vs tumor tissue in samples
#"NT" for normal tissue, "TP" for primoary tumor
library(TCGAbiolinks)
barcodes <- colnames(DataFilt)
TCGAquery_SampleTypes(barcodes, c("NB","NT","NBC","NEBV","NBM",))
TCGAquery_SampleTypes(barcodes, c("TP"))
clin <- GDCquery_clinic("TCGA-LUAD", type = "clinical", save.csv = TRUE)

gender_paiant=summarise(group_by(clin, gender),n = n())
barplot(muta)
genderbar<-ggplot(gender_paiant, aes(x=gender, y=n,color=gender)) +
  geom_bar(stat="identity",fill="white")+
  labs(title="Distribution of gender in TCGA",
       x="gender", 
       y = "number of gender in painten size")
genderbar+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))

race_paints=summarise(group_by(clin, race),n = n())
racebar<-ggplot(race_paints, aes(x=race, y=n,color=race)) +
  geom_bar(stat="identity",fill="white")+
  labs(title="Distribution of race in TCGA",
       x="race", 
       y = "number of race in painten size")
racebar+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))

vital_status_paints=summarise(group_by(clin, vital_status),n = n())
vital_statusbar<-ggplot(vital_status_paints, aes(x=vital_status, y=n,color=vital_status)) +
  geom_bar(stat="identity",fill="white")+
  labs(title="Distribution of vital_status in TCGA",
       x="vital_status", 
       y = "number of vital_status in painten size")
vital_statusbar+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))


age_at_diagnosis_paints=summarise(group_by(clin, age_at_diagnosis),n = n())
age_at_diagnosisbar<-ggplot(age_at_diagnosis_paints, aes(x=age_at_diagnosis, y=n,color=age_at_diagnosis)) +
  geom_bar(stat="identity",fill="white")+
  labs(title="Distribution of age_at_diagnosis in TCGA",
       x="age_at_diagnosis", 
       y = "number of age_at_diagnosis in painten size")
age_at_diagnosisbar+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))


ajcc_pathologic_tumor_stage_paints=summarise(group_by(clin, ajcc_pathologic_stage),n = n())
ajcc_pathologic_tumor_stagebar<-ggplot(ajcc_pathologic_tumor_stage_paints, aes(x=ajcc_pathologic_stage, y=n,color=ajcc_pathologic_stage)) +
  geom_bar(stat="identity",fill="white")+
  labs(title="Distribution of tumor stage in TCGA",
       x="tumor stage", 
       y = "number of tumor stage in painten size")
ajcc_pathologic_tumor_stagebar+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))



