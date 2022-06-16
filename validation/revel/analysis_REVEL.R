


library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggvenn)



#read files

lifted_DEG_MAF=read_csv("../input/lifted_DEG_MAF.csv")
drivers_moonlight=read_csv("../input/moonlight_driver.csv")
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


#moonlight

uniq_revel_moonlight_MAF=inner_join(drivers_moonlight,uniq_revel_DEG_MAF)

#anti subset
uniq_revel_anti_moonlight_MAF=anti_join(uniq_revel_DEG_MAF,uniq_revel_moonlight_MAF)


#means
mean_DEG=mean(uniq_revel_DEG_MAF$REVEL,na.rm=TRUE)
mean_moonlight=mean(uniq_revel_moonlight_MAF$REVEL,na.rm=TRUE)
mean_anti_moonlight=mean(uniq_revel_anti_moonlight_MAF$REVEL,na.rm=TRUE)



#ggplot histogram
df_moonlight_scores_analyse=data.frame(Name=rep("Moonlight",length(which(!is.na(uniq_revel_moonlight_MAF$REVEL)))),
                                       scores=uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$REVEL))])
df_DEG_scores_analyse=data.frame(Name=rep("DEG",length(which(!is.na(uniq_revel_DEG_MAF$REVEL)))),
                                 scores=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$REVEL))])
df_anti_moonlight_scores_analyse=data.frame(Name=rep("Anti_moonlight",length(which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL)))),
                                            scores=uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL))])

df_all=rbind(df_moonlight_scores_analyse,df_anti_moonlight_scores_analyse,df_DEG_scores_analyse)


mu=data.frame(Name=c("Moonlight","DEG","Anti_moonlight"),grp.mean=c(mean_moonlight,mean_DEG,mean_anti_moonlight))
ggplot_scores_all=ggplot(df_all, aes(x=scores,color=Name,fill=Name)) + 
  geom_histogram(binwidth=0.01, position = "identity", alpha = 0.32) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Name),linetype="dashed")+
  theme(legend.position="top")+
  labs(title="REVEL scores in coding region for\n three set with mean lines",
       x="REVEL scores for mutations in coding region", 
       y = "Count of mutations in scores bin")
ggplot_scores_all+theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))


#t.test for REVEL scores 
#revel scores for coding reagion
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))]))]
revel_score_uniq_moonlight=uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$Coding))][which(!is.na(uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$Coding))]))]
revel_score_uniq_anti_moonlight=uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$Coding))][which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$Coding))]))]

#running t.test
uniq_moonlight_vs_DEG=t.test(revel_score_uniq_moonlight,revel_score_uniq_DEG)
uniq_anti_moonlight_vs_DEG=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)
uniq_moonlight_vs_anti_moonlight=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)


print(uniq_moonlight_vs_DEG)
print(uniq_anti_moonlight_vs_DEG)
print(uniq_moonlight_vs_anti_moonlight)

#t.test for REVEL of all scores inlcuding non coding Cscape
#revel scores for coding reagion
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$REVEL))]
revel_score_uniq_moonlight=uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$REVEL))]
revel_score_uniq_anti_moonlight=uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL))]

#running t.test
uniq_moonlight_vs_DEG=t.test(revel_score_uniq_moonlight,revel_score_uniq_DEG)
uniq_anti_moonlight_vs_DEG=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)
uniq_moonlight_vs_anti_moonlight=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)

print(uniq_moonlight_vs_DEG)
print(uniq_anti_moonlight_vs_DEG)
print(uniq_moonlight_vs_anti_moonlight)


#revel scores for Coding reagion with threshold
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding>0.5))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding>0.5))]))]
revel_score_uniq_moonlight=uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$Coding>0.5))][which(!is.na(uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$Coding>0.5))]))]
revel_score_uniq_anti_moonlight=uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$Coding>0.5))][which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$Coding>0.5))]))]

#running t.test
uniq_moonlight_vs_DEG=t.test(revel_score_uniq_moonlight,revel_score_uniq_DEG)
uniq_anti_moonlight_vs_DEG=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)
uniq_moonlight_vs_anti_moonlight=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)



print(uniq_moonlight_vs_DEG)
print(uniq_anti_moonlight_vs_DEG)
print(uniq_moonlight_vs_anti_moonlight)


print(mean(revel_score_uniq_DEG))
print(mean(revel_score_uniq_moonlight))
print(mean(revel_score_uniq_anti_moonlight))




#revel scores for Coding reagion with threshold
revel_score_uniq_DEG=uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding>0.8))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding>0.8))]))]
revel_score_uniq_moonlight=uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$Coding>0.8))][which(!is.na(uniq_revel_moonlight_MAF$REVEL[which(!is.na(uniq_revel_moonlight_MAF$Coding>0.8))]))]
revel_score_uniq_anti_moonlight=uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$Coding>0.8))][which(!is.na(uniq_revel_anti_moonlight_MAF$REVEL[which(!is.na(uniq_revel_anti_moonlight_MAF$Coding>0.8))]))]

#running t.test
uniq_moonlight_vs_DEG=t.test(revel_score_uniq_moonlight,revel_score_uniq_DEG)
uniq_anti_moonlight_vs_DEG=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)
uniq_moonlight_vs_anti_moonlight=t.test(revel_score_uniq_anti_moonlight,revel_score_uniq_DEG)



print(uniq_moonlight_vs_DEG)
print(uniq_anti_moonlight_vs_DEG)
print(uniq_moonlight_vs_anti_moonlight)

coding_scores=length(which(!is.na(uniq_revel_DEG_MAF$Coding)))
revel_scores=length(which(!is.na(uniq_revel_DEG_MAF$Coding)))
missed_by_revel=length(which(is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))])))
missed_by_cscape=length(which(is.na(uniq_revel_DEG_MAF$Coding[which(!is.na(uniq_revel_DEG_MAF$REVEL))])))

t_run=length(revel_score_uniq_DEG)

#skatter plot

uniq_revel_DEG_MAF[which(!is.na(uniq_revel_DEG_MAF$Coding>0.8))][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding>0.8))]))][,"Revel","Coding"]
uniq_revel_DEG_MAF[which(!is.na(uniq_revel_DEG_MAF$Coding>0.8)),][which(!is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding>0.8))])),]
uniq_revel_DEG_MAF[,c("REVEL","Coding")]

coding_scores=length(which(!is.na(uniq_revel_DEG_MAF$Coding)))
revel_scores=length(which(!is.na(uniq_revel_DEG_MAF$Coding)))
missed_by_revel=length(which(is.na(uniq_revel_DEG_MAF$REVEL[which(!is.na(uniq_revel_DEG_MAF$Coding))])))
missed_by_cscape=length(which(is.na(uniq_revel_DEG_MAF$Coding[which(!is.na(uniq_revel_DEG_MAF$REVEL))])))

coding_uniq_revel_DEG_MAF=uniq_revel_DEG_MAF[,c("REVEL","Coding")][which(!is.na(uniq_revel_DEG_MAF$Coding)),]
REVEL_coding_uniq_revel_DEG_MAF=coding_uniq_revel_DEG_MAF[which(!is.na(coding_uniq_revel_DEG_MAF$REVEL)),]

point_REVELvsCscape=ggplot(REVEL_coding_uniq_revel_DEG_MAF, aes(x=REVEL, y=Coding)) +
  geom_point(size=0.01)+
  theme(legend.position="top")+
  labs(title="Revel scores vs Cscape coding region scores ",
       x="Revel scores", 
       y = "Cscape coding region scores")
point_REVELvsCscape +theme_linedraw()+theme(text = element_text(size = 20,family = "Times New Roman"))


REVEL_DEG_mut_Pr_gene_count=summarise(group_by(uniq_revel_DEG_MAF[which(!is.na(uniq_revel_DEG_MAF$REVEL)),], Hugo_Symbol),n = n())


REVEL_DEG_gene_mean=uniq_revel_DEG_MAF %>%
  group_by(Hugo_Symbol) %>%
  summarise_at(vars(REVEL), list(name = mean))



REVEL_moonlight_mut_Pr_gene_count=summarise(group_by(uniq_revel_moonlight_MAF[which(!is.na(uniq_revel_moonlight_MAF$REVEL)),], Hugo_Symbol),n = n())


REVEL_moonlight_gene_mean=uniq_revel_moonlight_MAF %>%
  group_by(Hugo_Symbol) %>%
  summarise_at(vars(REVEL), list(name = mean))

