#REVEL analysys

#load libraties
library(tidyverse)
library(patchwork)
library(ggplot2)

#read files

lifted_DEG_MAF=read_csv("../input/lifted_DEG_MAF.csv")
drivers_moonlight=read_csv("../input/moonlight_driver.csv")
revel_DEG=read_csv("./scores_DEG")
results_DEG=read_tsv("../../Mut_CancerDriver/results/cscape_DEG.txt")

#mutate data, remove extra index and add chr to chromosome
revel_DEG=dplyr::select(revel_DEG,-1)
revel_DEG$X1=paste("chr",revel_DEG$X1,sep="")

results_DEG=dplyr::rename(results_DEG,Chromosome="# Chromosome")
results_DEG=dplyr::rename(results_DEG,Start_Position=Position)
results_DEG=dplyr::rename(results_DEG,Tumor_Seq_Allele2=Mutant)
results_DEG=dplyr::rename(results_DEG,Reference_Allele=Reference)
results_DEG$Chromosome=paste("chr",results_DEG$Chromosome,sep="")

#find uniq entries
uniq_lifted_DEG_MAF=distinct(lifted_DEG_MAF,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)
uniq_revel_DEG=distinct(revel_DEG,X1,X2,X3,X4, .keep_all = TRUE)
uniq_results_DEG=distinct(results_DEG,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)


#join maf to results
uniq_revel_DEG_MAF2=left_join(uniq_revel_DEG,uniq_lifted_DEG_MAF,by=c(X1="Chromosome",X2= "Start_Position",X3 ="Reference_Allele", X4="Tumor_Seq_Allele2"))
uniq_revel_DEG_MAF=left_join(uniq_revel_DEG_MAF2,uniq_results_DEG,by=c(X1="Chromosome",X2= "Start_Position",X3 ="Reference_Allele", X4="Tumor_Seq_Allele2"))


#subset for moonlight drivers
uniq_revel_moonlight_MAF=inner_join(drivers_moonlight,uniq_revel_DEG_MAF)

#anti subset
uniq_revel_anti_moonlight_MAF=anti_join(uniq_revel_DEG_MAF,uniq_revel_moonlight_MAF)



#ggplot histogram
#create data frame
df_moonlight_scores_analyse=data.frame(name=rep("Moonlight",length(which(!is.na(uniq_revel_moonlight_MAF$REVEL)))),
                                       scores=uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$REVEL))])
df_DEG_scores_analyse=data.frame(name=rep("DEG",length(which(!is.na(uniq_revel_DEG_MAF$REVEL)))),
                                 scores=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$REVEL))])
df_anti_moonlight_scores_analyse=data.frame(name=rep("Anti_moonlight",length(which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL)))),
                                            scores=uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL))])

df_all=rbind(df_moonlight_scores_analyse,df_anti_moonlight_scores_analyse,df_DEG_scores_analyse)

#means
mean_DEG=mean(uniq_revel_DEG_MAF$REVEL,na.rm=TRUE)
mean_moonlight=mean(uniq_revel_moonlight_MAF$REVEL,na.rm=TRUE)
mean_anti_moonlight=mean(uniq_revel_anti_moonlight_MAF$REVEL,na.rm=TRUE)
mu=data.frame(name=c("Moonlight","DEG","Anti_moonlight"),grp.mean=c(mean_moonlight,mean_DEG,mean_anti_moonlight))

#ploting histogram
ggplot_scores_all=ggplot(df_all, aes(x=scores,color=name)) + geom_histogram(fill="white",binwidth=0.01) +geom_vline(data=mu, aes(xintercept=grp.mean, color=name),linetype="dashed")+theme(legend.position="top")+labs(title="Cscape histrogram of coding scores",x="Scores coding region", y = "Count")
ggplot_scores_all

