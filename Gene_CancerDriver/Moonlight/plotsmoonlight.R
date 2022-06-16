#Moonlight plots


#library(UpSetR)

jpeg("BarplotMoonlightGenes2.png", width = 350, height = 350)
#loads list of found driver genes
DataTS <- read.csv("../results/Moon_CancerGenes_EB_TSG.csv")
DataOC <- read.csv("../results/Moon_CancerGenes_EB_OCG.csv")

#countshow many genese defined as both a tumoor supprosr and an onco gene by moonlight

Dualcount <- length(intersect(DataTS$x,DataOC$x))
TScount <- length(DataTS$x)
OCcount <- length(DataOC$x)



#creates an input list for Upset
#listinput <- list("Tomor Suppresor"=c(rep(1,length(DataTS$x)-dual),rep(3,dual)), "onco genes"=c(rep(2,length(DataOC$x)-dual),rep(3,dual)))
#plotmoon=upset(fromList(listinput))




GeneDriverBarPlot=barplot(c(TScount,OCcount),
        main = "Moonlight predictions",
        xlab = "Gene count",
        ylab = "Driver genes",
        names.arg = c("Tumor suppressors","Onco genes"),
        col = c("green","red"),
        horiz = TRUE)
#ggsave(filename = "BarplotMoonlightGenes.pdf")


dev.off()


