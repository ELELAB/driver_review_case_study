#Analay and format for result files

library(tidyverse)

#load data

lifted_moonlight_MAF=read_csv("./lifted_moonlight_MAF.csv")
cscape_input_moonlight_MAF=read_csv("./cscape_input_converted_MAF_moonlight.csv",col_names=c("Chromosome",))
results_moonlight=read_tsv("../results/cscape_mooonlight.txt")

results_moonlight=dplyr::rename(results_moonlight,Chromosome="# Chromosome")
results_moonlight=dplyr::rename(results_moonlight,Start_Position=Position)
results_moonlight=dplyr::rename(results_moonlight,Tumor_Seq_Allele2=Mutant)
results_moonlight=dplyr::rename(results_moonlight,Reference_Allele=Reference)


"modifly "
results_moonlight$Chromosome=paste("chr",results_moonlight$Chromosome,sep="")



results_moonlight_MAF=right_join(results_moonlight,lifted_moonlight_MAF,by=c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))

fiskMAF3=group_by(results_moonlight_MAF,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2)
fiskMAF4=summarise(fisk3)

lifted_DEG_MAF=read_csv("./lifted_DEG_MAF.csv")
cscape_input_DEG_MAF=read_csv("./cscape_input_converted_MAF_DEG.csv")
results_DEG=read_tsv("../results/cscape_DEG.txt")

results_DEG=dplyr::rename(results_DEG,Chromosome="# Chromosome")
results_DEG=dplyr::rename(results_DEG,Start_Position=Position)


results_DEG$Chromosome=paste("chr",results_DEG$Chromosome,sep="")


results_DEG_MAF=left_join(results_DEG,lifted_DEG_MAF)

hist(results_DEG$Coding)
hist(results_moonlight$Coding)
x1=mean(results_DEG$Coding,na.rm=TRUE)
x2=mean(results_moonlight$Coding,na.rm=TRUE)

t.test(x1,x2)
x3=t.test(results_DEG$Coding,results_moonlight$Coding)

cscape_driver_moonlight_total=length(which(!is.na(results_moonlight$Coding)==TRUE))

cscape_driver_moonlight_05=length(which(results_moonlight$Coding>0.5))

cscape_driver_moonlight_08=length(which(results_moonlight$Coding>0.8))

frac_moonlight_05=cscape_driver_moonlight_05/cscape_driver_moonlight_total
frac_moonlight_08=cscape_driver_moonlight_08/cscape_driver_moonlight_total

uniq_gene_moonlight=



cscape_driver_DEG_total=length(which(!is.na(results_DEG$Coding)==TRUE))

cscape_driver_DEG_05=length(which(results_DEG$Coding>0.5))

cscape_driver_DEG_08=length(which(results_DEG$Coding>0.8))

frac_DEG_05=cscape_driver_DEG_05/cscape_driver_DEG_total
frac_DEG_08=cscape_driver_DEG_08/cscape_driver_DEG_total

fisk=anti_join(results_moonlight_MAF,results_moonlight,by=)

anti_DEG_driverMut=anti_join(results_DEG_MAF,results_moonlight_MAF)
anti_moonlight_driverMut=anti_join(results_DEG_MAF,results_moonlight_MAF)

chromosome Start_position
