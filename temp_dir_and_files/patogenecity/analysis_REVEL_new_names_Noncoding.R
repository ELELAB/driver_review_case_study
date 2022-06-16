


library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggvenn)



#read files

lifted_DEG_MAF=read_csv("../input/lifted_DEG_MAF.csv")
drivers_Oncogenic_mediator=read_csv("../input/moonlight_driver.csv")
revel_DEG=read_csv("./scores_DEG")


results_DEG=read_tsv("../../Mut_CancerDriver/results/cscape_DEG.txt")

#remove extra colum
revel_DEG=dplyr::select(revel_DEG,-1)
revel_DEG$X1=paste("chr",revel_DEG$X1,sep="")

results_DEG=dplyr::rename(results_DEG,Chromosome="# Chromosome")
results_DEG=dplyr::rename(results_DEG,Start_Position=Position)
results_DEG=dplyr::rename(results_DEG,Tumor_Seq_Allele2=Mutant)
results_DEG=dplyr::rename(results_DEG,Reference_Allele=Reference)
results_DEG$Chromosome=paste("chr",results_DEG$Chromosome,sep="")
#find uniq entryes
uniq_lifted_DEG_MAF=distinct(lifted_DEG_MAF,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)
uniq_revel_DEG=distinct(revel_DEG,X1,X2,X3,X4, .keep_all = TRUE)
uniq_results_DEG=distinct(results_DEG,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2, .keep_all = TRUE)

#join MAF data

uniq_revel_DEG_MAF2=left_join(uniq_revel_DEG,uniq_lifted_DEG_MAF,by=c(X1="Chromosome",X2= "Start_Position",X3 ="Reference_Allele", X4="Tumor_Seq_Allele2"))
uniq_revel_DEG_MAF=left_join(uniq_revel_DEG_MAF2,uniq_results_DEG,by=c(X1="Chromosome",X2= "Start_Position",X3 ="Reference_Allele", X4="Tumor_Seq_Allele2"))


#Oncogenic_mediator

uniq_revel_Oncogenic_mediator_MAF=inner_join(drivers_Oncogenic_mediator,uniq_revel_DEG_MAF)

#anti subset
uniq_revel_None_Driver_MAF=anti_join(uniq_revel_DEG_MAF,uniq_revel_Oncogenic_mediator_MAF)


#means
mean_DEG=mean(uniq_revel_DEG_MAF$REVEL,na.rm=TRUE)
mean_Oncogenic_mediator=mean(uniq_revel_Oncogenic_mediator_MAF$REVEL,na.rm=TRUE)
mean_None_Driver=mean(uniq_revel_None_Driver_MAF$REVEL,na.rm=TRUE)



#ggplot histogram
df_Oncogenic_mediator_scores_analyse=data.frame(Name=rep("Oncogenic mediator",length(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))]))])),
                                                 scores=uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))]))])
df_DEG_scores_analyse=data.frame(Name=rep("DEG",length(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))]))])),
                                 scores=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))]))])
df_None_Driver_scores_analyse=data.frame(Name=rep("None Driver",length(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))]))])),
                                         scores=uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))]))])


df_all=rbind(df_Oncogenic_mediator_scores_analyse,df_None_Driver_scores_analyse,df_DEG_scores_analyse)


mu=data.frame(Name=c("Oncogenic mediator","DEG","None Driver"),grp.mean=c(mean_Oncogenic_mediator,mean_DEG,mean_None_Driver))
ggplot_scores_all=ggplot(df_all, aes(x=scores,color=Name,fill=Name)) + 
  geom_histogram(binwidth=0.05, position = "identity", alpha = 0.32) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Name),linetype="dashed")+
  theme(legend.position="top")+
  labs(title="REVEL scores in Noncoding region for\n three set with mean lines",
       x="REVEL scores for mutations in Noncoding region", 
       y = "Count of mutations in scores bin")
ggplot_scores_all+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))

#ggplot histogram log noncoding
df_Oncogenic_mediator_scores_analyse=data.frame(Name=rep("Oncogenic mediator",length(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))]))])),
                                                scores=log(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))]))]))
df_DEG_scores_analyse=data.frame(Name=rep("DEG",length(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))]))])),
                                 scores=log(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))]))]))
df_None_Driver_scores_analyse=data.frame(Name=rep("None Driver",length(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))]))])),
                                         scores=log(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))]))]))

df_all=rbind(df_Oncogenic_mediator_scores_analyse,df_None_Driver_scores_analyse,df_DEG_scores_analyse)
mean_DEG=mean(log(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))]))]))
mean_Oncogenic_mediator=mean(log(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))]))]))
mean_None_Driver=mean(log(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))]))]))


mu=data.frame(Name=c("Oncogenic mediator","DEG","None Driver"),grp.mean=c(mean_Oncogenic_mediator,mean_DEG,mean_None_Driver))
ggplot_scores_all=ggplot(df_all, aes(x=scores,color=Name,fill=Name)) + 
  geom_histogram(binwidth=0.4, position = "identity", alpha = 0.32) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Name),linetype="dashed")+
  theme(legend.position="top")+
  labs(title="REVEL log transformed scores in Noncoding region for\n three set with mean lines",
       x="REVEL scores for mutations in Noncoding region", 
       y = "Count of mutations in scores bin")
ggplot_scores_all+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))



#ggplot histogram coding
df_Oncogenic_mediator_scores_analyse=data.frame(Name=rep("Oncogenic mediator",length(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Coding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Coding))]))])),
                                                scores=uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Coding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Coding))]))])
df_DEG_scores_analyse=data.frame(Name=rep("DEG",length(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))]))])),
                                 scores=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))]))])
df_None_Driver_scores_analyse=data.frame(Name=rep("None Driver",length(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Coding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Coding))]))])),
                                         scores=uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Coding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Coding))]))])

df_all=rbind(df_Oncogenic_mediator_scores_analyse,df_None_Driver_scores_analyse,df_DEG_scores_analyse)


mu=data.frame(Name=c("Oncogenic mediator","DEG","None Driver"),grp.mean=c(mean_Oncogenic_mediator,mean_DEG,mean_None_Driver))
ggplot_scores_all=ggplot(df_all, aes(x=scores,color=Name,fill=Name)) + 
  geom_histogram(binwidth=0.01, position = "identity", alpha = 0.32) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Name),linetype="dashed")+
  theme(legend.position="top")+
  labs(title="REVEL scores in coding region for\n three set with mean lines",
       x="REVEL scores for mutations in coding region", 
       y = "Count of mutations in scores bin")
ggplot_scores_all+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))

#t.test for REVEL scores 
#revel scores for Noncoding reagion
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))]))]
revel_score_uniq_Oncogenic_mediator=uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding))]))]
revel_score_uniq_None_Driver=uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding))]))]

#running t.test
uniq_Oncogenic_mediator_vs_DEG=t.test(revel_score_uniq_Oncogenic_mediator,revel_score_uniq_DEG)
uniq_None_Driver_vs_DEG=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)
uniq_Oncogenic_mediator_vs_None_Driver=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)


print(uniq_Oncogenic_mediator_vs_DEG)
print(uniq_None_Driver_vs_DEG)
print(uniq_Oncogenic_mediator_vs_None_Driver)

#t.test for REVEL of all scores inlcuding non Noncoding Cscape
#revel scores for Noncoding reagion
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$REVEL))]
revel_score_uniq_Oncogenic_mediator=uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL))]
revel_score_uniq_None_Driver=uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$REVEL))]

#running t.test
uniq_Oncogenic_mediator_vs_DEG=t.test(revel_score_uniq_Oncogenic_mediator,revel_score_uniq_DEG)
uniq_None_Driver_vs_DEG=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)
uniq_Oncogenic_mediator_vs_None_Driver=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)

print(uniq_Oncogenic_mediator_vs_DEG)
print(uniq_None_Driver_vs_DEG)
print(uniq_Oncogenic_mediator_vs_None_Driver)


#revel scores for Noncoding reagion with threshold
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.5))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.5))]))]
revel_score_uniq_Oncogenic_mediator=uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding>0.5))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding>0.5))]))]
revel_score_uniq_None_Driver=uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding>0.5))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding>0.5))]))]

#running t.test
uniq_Oncogenic_mediator_vs_DEG=t.test(revel_score_uniq_Oncogenic_mediator,revel_score_uniq_DEG)
uniq_None_Driver_vs_DEG=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)
uniq_Oncogenic_mediator_vs_None_Driver=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)



print(uniq_Oncogenic_mediator_vs_DEG)
print(uniq_None_Driver_vs_DEG)
print(uniq_Oncogenic_mediator_vs_None_Driver)


print(mean(revel_score_uniq_DEG))
print(mean(revel_score_uniq_Oncogenic_mediator))
print(mean(revel_score_uniq_None_Driver))




#revel scores for Noncoding reagion with threshold
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.8))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.8))]))]
revel_score_uniq_Oncogenic_mediator=uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding>0.8))][which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$Noncoding>0.8))]))]
revel_score_uniq_None_Driver=uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding>0.8))][which(!is.na(uniq_revel_None_Driver_MAF$REVEL[which(!is.na(uniq_revel_None_Driver_MAF$Noncoding>0.8))]))]

#running t.test
uniq_Oncogenic_mediator_vs_DEG=t.test(revel_score_uniq_Oncogenic_mediator,revel_score_uniq_DEG)
uniq_None_Driver_vs_DEG=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)
uniq_Oncogenic_mediator_vs_None_Driver=t.test(revel_score_uniq_None_Driver,revel_score_uniq_DEG)



print(uniq_Oncogenic_mediator_vs_DEG)
print(uniq_None_Driver_vs_DEG)
print(uniq_Oncogenic_mediator_vs_None_Driver)

Noncoding_scores=length(which(!is.na(uniq_revel_DEG_MAF$Noncoding)))
revel_scores=length(which(!is.na(uniq_revel_DEG_MAF$REVEL)))
missed_by_revel=length(which(is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))])))
missed_by_cscape=length(which(is.na(uniq_revel_DEG_MAF$Noncoding[which(!is.na(uniq_revel_DEG_MAF$REVEL))])))

t_run=length(revel_score_uniq_DEG)

#skatter plot

uniq_revel_DEG_MAF[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.8))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.8))]))][,"Revel","Noncoding"]
uniq_revel_DEG_MAF[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.8)),][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding>0.8))])),]
uniq_revel_DEG_MAF[,c("REVEL","Noncoding")]

Noncoding_scores=length(which(!is.na(uniq_revel_DEG_MAF$Noncoding)))
revel_scores=length(which(!is.na(uniq_revel_DEG_MAF$Noncoding)))
missed_by_revel=length(which(is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Noncoding))])))
missed_by_cscape=length(which(is.na(uniq_revel_DEG_MAF$Noncoding[which(!is.na(uniq_revel_DEG_MAF$REVEL))])))

Noncoding_uniq_revel_DEG_MAF=uniq_revel_DEG_MAF[,c("REVEL","Noncoding")][which(!is.na(uniq_revel_DEG_MAF$Noncoding)),]
REVEL_Noncoding_uniq_revel_DEG_MAF=Noncoding_uniq_revel_DEG_MAF[which(!is.na(Noncoding_uniq_revel_DEG_MAF$REVEL)),]

point_REVELvsCscape=ggplot(REVEL_Noncoding_uniq_revel_DEG_MAF, aes(x=REVEL, y=Noncoding)) +
  geom_point(size=0.01)+
  theme(legend.position="top")+
  labs(title="Revel scores vs Cscape Noncoding region scores ",
       x="Revel scores", 
       y = "Cscape Noncoding region scores")
point_REVELvsCscape +theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))


REVEL_DEG_mut_Pr_gene_count=summarise(group_by(uniq_revel_DEG_MAF[which(!is.na(uniq_revel_DEG_MAF$REVEL)),], Hugo_Symbol),n = n())


REVEL_DEG_gene_mean=uniq_revel_DEG_MAF %>%
  group_by(Hugo_Symbol) %>%
  summarise_at(vars(REVEL), list(name = mean))



REVEL_Oncogenic_mediator_mut_Pr_gene_count=summarise(group_by(uniq_revel_Oncogenic_mediator_MAF[which(!is.na(uniq_revel_Oncogenic_mediator_MAF$REVEL)),], Hugo_Symbol),n = n())


REVEL_Oncogenic_mediator_gene_mean=uniq_revel_Oncogenic_mediator_MAF %>%
  group_by(Hugo_Symbol) %>%
  summarise_at(vars(REVEL), list(name = mean))

