##R script for MAF convertions to Cscape input

#load functions
library(tidyverse)
source("/data/user/shared_projects/moonlight2/mutations/old_src/small_functions.R")

#Read data


input_MAF=read_csv("../inputfiles/mutations.csv")
DEG=get(load("../../Gene_CancerDriver/inputfiles/LUAD_dataDEGs_deconvolution.rda"))


#removes extra colum and adds driver colum 
DEG_gene=row.names(DEG)
Frame_DEG_gene=tibble(DEG_gene)
Frame_DEG_gene=dplyr::rename(Frame_DEG_gene,SYMBOL=DEG_gene)

write_csv(Frame_DEG_gene,file="./Frame_DEG_gene.csv")

#MAF with genes from moonlight
DEG_MAF=inner_join(Frame_DEG_gene,input_MAF)
write_csv(DEG_MAF,file="./DEG_MAF.csv")


#lift GRCh38
lifted_MAF=LiftMAF(DEG_MAF,"GRCh38")
write_csv(lifted_MAF,file="./lifted_DEG_MAF.csv")

#convert MAF to cscape input

cscape_input=MAFtoCscape(lifted_MAF)

#rite file
write.table(cscape_input,"./cscape_input_converted_MAF_DEG.csv",sep=",",col.names=FALSE,row.names=FALSE,quote = FALSE)

#write_csv(cscape_input,file="./cscape_input_converted_MAF_DEG.csv")


