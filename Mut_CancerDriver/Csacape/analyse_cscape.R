#Analay and format for result files

library(tidyverse)
library(patchwork)
library(ggplot2)

#load data

lifted_moonlight_MAF=read_csv("./lifted_moonlight_MAF.csv")
cscape_input_moonlight_MAF=read_csv("./cscape_input_converted_MAF_moonlight.csv",col_names=c("Chromosome",))
results_moonlight=read_tsv("../results/cscape_mooonlight.txt")


#modify moonlug tresult file to 
results_moonlight=dplyr::rename(results_moonlight,Chromosome="# Chromosome")
results_moonlight=dplyr::rename(results_moonlight,Start_Position=Position)
results_moonlight=dplyr::rename(results_moonlight,Tumor_Seq_Allele2=Mutant)
results_moonlight=dplyr::rename(results_moonlight,Reference_Allele=Reference)
results_moonlight$Chromosome=paste("chr",results_moonlight$Chromosome,sep="")

#find unice entries
uniq_results_moonlight=distinct(results_moonlight,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)
uniq_lifted_moonlight_MAF=distinct(lifted_moonlight_MAF,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)

#join mooonligh result with MAF file
uniq_results_moonlight_MAF=left_join(uniq_results_moonlight,uniq_lifted_moonlight_MAF,by=c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))

#find doblications
doblicates_results_moonlight=results_moonlight[duplicated(results_moonlight),]



#load data
lifted_DEG_MAF=read_csv("./lifted_DEG_MAF.csv")
cscape_input_DEG_MAF=read_csv("./cscape_input_converted_MAF_DEG.csv")
results_DEG=read_tsv("../results/cscape_DEG.txt")

#modify results
results_DEG=dplyr::rename(results_DEG,Chromosome="# Chromosome")
results_DEG=dplyr::rename(results_DEG,Start_Position=Position)
results_DEG=dplyr::rename(results_DEG,Tumor_Seq_Allele2=Mutant)
results_DEG=dplyr::rename(results_DEG,Reference_Allele=Reference)
results_DEG$Chromosome=paste("chr",results_DEG$Chromosome,sep="")

#find unice entries
uniq_results_DEG=distinct(results_DEG,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)
uniq_lifted_DEG_MAF=distinct(lifted_DEG_MAF,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)

#join DEG result with MAF file
uniq_results_DEG_MAF=left_join(uniq_results_DEG,uniq_lifted_DEG_MAF,by=c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))

#find doblications
doblicates_results_DEG=results_DEG[duplicated(results_DEG),]



#anti join (anti moonlight)
uniq_results_anti_moonlight_MAF=anti_join(uniq_results_DEG_MAF,uniq_results_moonlight_MAF)


###histogram
#means
mean_DEG=mean(results_DEG$Coding,na.rm=TRUE)
mean_moonlight=mean(results_moonlight$Coding,na.rm=TRUE)
mean_anti_moonlight=mean(uniq_results_anti_moonlight_MAF$Coding,na.rm=TRUE)


#ggplot histogram
df_moonlight_scores_analyse=data.frame(Name=rep("Oncogenic mediators",length(which(!is.na(uniq_results_moonlight_MAF$Coding)))),
                                       scores=uniq_results_moonlight_MAF$Coding[which(!is.na(uniq_results_moonlight_MAF$Coding))])
df_DEG_scores_analyse=data.frame(Name=rep("DEG",length(which(!is.na(uniq_results_DEG_MAF$Coding)))),
                                 scores=uniq_results_DEG_MAF$Coding[which(!is.na(uniq_results_DEG_MAF$Coding))])
df_anti_moonlight_scores_analyse=data.frame(Name=rep("Non-Drivers",length(which(!is.na(uniq_results_anti_moonlight_MAF$Coding)))),
                                            scores=uniq_results_anti_moonlight_MAF$Coding[which(!is.na(uniq_results_anti_moonlight_MAF$Coding))])

df_all=rbind(df_moonlight_scores_analyse,df_anti_moonlight_scores_analyse,df_DEG_scores_analyse)


#ggplot all
#mu <- ddply(df_all, "name", summarise, grp.mean=mean(name))
mu=data.frame(Name=c("Oncogenic mediators","Non-Drivers","DEG"),grp.mean=c(mean_moonlight,mean_anti_moonlight,mean_DEG))
ggplot_scores_all=ggplot(df_all, aes(x=scores,fill=Name,color=Name)) + 
  geom_histogram(binwidth=0.01, position = "identity", alpha = 0.12) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Name),linetype="dashed")+
  theme(legend.position="top")+
  labs(title="Cscape scores in coding region for\n three set with mean lines",
       x="Cscape scores for mutations in coding region", 
       y = "count of mutations in scores bin")
ggplot_scores_all+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))


###t.test
uniq_moonlight_vs_DEG=t.test(uniq_results_DEG_MAF$Coding,uniq_results_moonlight_MAF$Coding)
uniq_anti_moonlight_vs_DEG=t.test(uniq_results_DEG_MAF$Coding,uniq_results_anti_moonlight_MAF$Coding)
uniq_moonlight_vs_anti_moonlight=t.test(uniq_results_moonlight_MAF$Coding,uniq_results_anti_moonlight_MAF$Coding)


###frequency of hits for resutls moonlinglit in coding region
#moonlight
cscape_driver_mut_moonlight_total=length(which(!is.na(uniq_results_moonlight_MAF$Coding)==TRUE))
cscape_driver_mut_moonlight_05=length(which(uniq_results_moonlight_MAF$Coding>0.5))
cscape_driver_mut_moonlight_08=length(which(uniq_results_moonlight_MAF$Coding>0.8))

freq_moonlight_05=cscape_driver_mut_moonlight_05/cscape_driver_mut_moonlight_total
freq_moonlight_08=cscape_driver_mut_moonlight_08/cscape_driver_mut_moonlight_total

#DEG
cscape_driver_mut_DEG_total=length(which(!is.na(uniq_results_DEG_MAF$Coding)==TRUE))
cscape_driver_mut_DEG_05=length(which(uniq_results_DEG_MAF$Coding>0.5))
cscape_driver_mut_DEG_08=length(which(uniq_results_DEG_MAF$Coding>0.8))

freq_DEG_05=cscape_driver_mut_DEG_05/cscape_driver_mut_DEG_total
freq_DEG_08=cscape_driver_mut_DEG_08/cscape_driver_mut_DEG_total

#anti moonlight
cscape_driver_mut_anti_moonlight_total=length(which(!is.na(uniq_results_anti_moonlight_MAF$Coding)==TRUE))
cscape_driver_mut_anti_moonlight_05=length(which(uniq_results_anti_moonlight_MAF$Coding>0.5))
cscape_driver_mut_anti_moonlight_08=length(which(uniq_results_anti_moonlight_MAF$Coding>0.8))

freq_anti_moonlight_05=cscape_driver_mut_anti_moonlight_05/cscape_driver_mut_anti_moonlight_total
freq_anti_moonlight_08=cscape_driver_mut_anti_moonlight_08/cscape_driver_mut_anti_moonlight_total



###Fraction hits in coding region per gene
#moonlight
#uniq gene
cscape_driver_mut_DEG_total=length(which(!is.na(uniq_results_DEG_MAF$Coding)==TRUE))
gene_results_moonlight=uniq_results_moonlight_MAF$Hugo_Symbol[which(!is.na(uniq_results_moonlight_MAF$Coding))]
uniq_gene_results_moonlight=unique(gene_results_moonlight)

#hits over threshold per gene in threshold
fraction_moonlight=cscape_driver_mut_moonlight_total/length(uniq_gene_results_moonlight)

#hits over treshhold to all genes
fraction_moonlight_thres05_vsALL=cscape_driver_mut_moonlight_05/length(uniq_gene_results_moonlight)
fraction_moonlight_thres08_vsALL=cscape_driver_mut_moonlight_08/length(uniq_gene_results_moonlight)


#DEG
#uniq gene 
gene_results_DEG=uniq_results_DEG_MAF$Hugo_Symbol[which(!is.na(uniq_results_DEG_MAF$Coding))]
uniq_gene_results_DEG=unique(gene_results_DEG)

#hits over threshold per gene 
fraction_DEG=cscape_driver_mut_DEG_total/length(uniq_gene_results_DEG)

#hits over treshhold to all genes
fraction_DEG_thres05_vsALL=cscape_driver_mut_DEG_05/length(uniq_gene_results_DEG)
fraction_DEG_thres08_vsALL=cscape_driver_mut_DEG_08/length(uniq_gene_results_DEG)

#!#case without KRAS

cscape_driver_mut_DEG_total_minusKRAS=length(which(!is.na(uniq_results_DEG_MAF[which(uniq_results_DEG_MAF$Hugo_Symbol!="KRAS"),]$Coding)==TRUE))
cscape_driver_mut_DEG_05_minusKRAS=length(which(uniq_results_DEG_MAF[which(uniq_results_DEG_MAF$Hugo_Symbol!="KRAS"),]$Coding>0.5))
cscape_driver_mut_DEG_08_minusKRAS=length(which(uniq_results_DEG_MAF[which(uniq_results_DEG_MAF$Hugo_Symbol!="KRAS"),]$Coding>0.8))

fraction_DEG_minusKRAS=cscape_driver_mut_DEG_total_minusKRAS/(length(uniq_gene_results_DEG)-1)
fraction_DEG_thres05_minusKRAS=cscape_driver_mut_DEG_05_minusKRAS/(length(uniq_gene_results_DEG)-1)
fraction_DEG_thres08_minusKRAS=cscape_driver_mut_DEG_08_minusKRAS/(length(uniq_gene_results_DEG)-1)



#anti_moonlight
#uniq gene in theshold
gene_results_anti_moonlight=uniq_results_anti_moonlight_MAF$Hugo_Symbol[which(!is.na(uniq_results_anti_moonlight_MAF$Coding))]
uniq_gene_results_anti_moonlight=unique(gene_results_anti_moonlight)

#hits over threshold
fraction_anti_moonlight=cscape_driver_mut_anti_moonlight_total/length(uniq_gene_results_anti_moonlight)

#hits over treshhold to all genes
fraction_anti_moonlight_thres05_vsALL=cscape_driver_mut_anti_moonlight_05/length(uniq_gene_results_anti_moonlight)
fraction_anti_moonlight_thres08_vsALL=cscape_driver_mut_anti_moonlight_08/length(uniq_gene_results_anti_moonlight)



###hist of mutation pr gene distribution

#moonlight
fisk=uniq_results_moonlight_MAF[which(!is.na(uniq_results_moonlight_MAF$Coding)),]
hits_per_gene_results_moonlight=summarise(group_by(uniq_results_moonlight_MAF[which(!is.na(uniq_results_moonlight_MAF$Coding)),], Hugo_Symbol),n = n())
hist_dist_mut_driver_moonlight=ggplot(hits_per_gene_results_moonlight, aes(x=n)) + 
  geom_histogram(binwidth=1)+
  theme(legend.position="top")+
  labs(title="Moonlight mutation\n in coding region pr gene",
       x="Number of mutation in a gene", 
       y = "Number of genes with\n mutations in coding region")+
  xlim(0,25)


#DEG
hits_per_gene_results_DEG=summarise(group_by(uniq_results_DEG_MAF[which(!is.na(uniq_results_DEG_MAF$Coding)),], Hugo_Symbol),n = n())
hist_dist_mut_driver_DEG=ggplot(hits_per_gene_results_DEG, aes(x=n)) +
  geom_histogram(binwidth=1)+
  theme(legend.position="top")+
  labs(title="DEG mutation\n in coding region pr gene",
       x="Number of mutation in a gene", 
       y = "Number of genes with\n mutations in coding region")+
  xlim(0,25)


#DEG without KRAS
uniq_results_DEG_MAF_without_KRAS=uniq_results_DEG_MAF[which(uniq_results_DEG_MAF$Hugo_Symbol!="KRAS"),]
hits_per_gene_results_DEG=summarise(group_by(uniq_results_DEG_MAF_without_KRAS[which(!is.na(uniq_results_DEG_MAF_without_KRAS$Coding)),], Hugo_Symbol),n = n())
hist_dist_mut_driver_DEG_minus_KRAS=ggplot(hits_per_gene_results_DEG, aes(x=n)) +
  geom_histogram(binwidth=1)+
  theme(legend.position="top")+
  labs(title="DEG without KRAS mutation\n in coding region pr gene",
       x="Number of mutation in a gene", 
       y = "Number of genes with\n mutations in coding region")

#anti_moonlight
hits_per_gene_results_anti_moonlight=summarise(group_by(uniq_results_anti_moonlight_MAF[which(!is.na(uniq_results_anti_moonlight_MAF$Coding)),], Hugo_Symbol),n = n())
hist_dist_mut_driver_anti_moonlight=ggplot(hits_per_gene_results_anti_moonlight, aes(x=n)) + 
  geom_histogram(binwidth=1)+
  theme(legend.position="top")+
  labs(title="Anti-Moonlight mutation in\n coding region pr gene",
       x="Number of mutation in a gene", 
       y = "Number of genes with\n mutations in coding region")+
  xlim(0,25)


#anti_moonlight without KRAS
uniq_results_anti_moonlight_MAF_without_KRAS=uniq_results_anti_moonlight_MAF[which(uniq_results_anti_moonlight_MAF$Hugo_Symbol!="KRAS"),]
hits_per_gene_results_anti_moonlight=summarise(group_by(uniq_results_anti_moonlight_MAF_without_KRAS[which(!is.na(uniq_results_anti_moonlight_MAF_without_KRAS$Coding)),], Hugo_Symbol),n = n())
hist_dist_mut_driver_anti_moonlight_mius_KRAS=ggplot(hits_per_gene_results_anti_moonlight, aes(x=n)) +
  geom_histogram(binwidth=1)+
  theme(legend.position="top")+
  labs(title="Anti-Moonlight without KRAS\n mutation in coding region pr gene",
       x="Number of mutation in a gene", 
       y = "Number of genes with\n mutations in coding region")


hist_dist_mut_driver_moonlight/hist_dist_mut_driver_DEG/hist_dist_mut_driver_anti_moonlight


hist_dist_mut_driver_moonlight/hist_dist_mut_driver_DEG+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))/hist_dist_mut_driver_anti_moonlight+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))
(hist_dist_mut_driver_moonlight+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))/
  hist_dist_mut_driver_DEG+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))/
  hist_dist_mut_driver_anti_moonlight+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman")))


#mutations type

mutations=summarise(group_by(uniq_lifted_DEG_MAF, Variant_Type),n = n())
barplot(muta)
mutbar<-ggplot(mutations, aes(x=Variant_Type, y=n,color=Variant_Type)) +
  geom_bar(stat="identity",fill="white")+
  labs(title="Distribution of mutation types in MAF file",
       x="mutations type", 
       y = "Number of mutation by type")
mutbar+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))
