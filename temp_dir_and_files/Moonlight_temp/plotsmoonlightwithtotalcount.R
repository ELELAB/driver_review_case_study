#Moonlight plots


#library(UpSetR)

jpeg("BarplotMoonlightGenes3.png", width = 350, height = 350)
#loads list of found driver genes
DataTS <- read.csv("../results/Moon_CancerGenes_EB_TSG.csv")
DataOC <- read.csv("../results/Moon_CancerGenes_EB_OCG.csv")
dataDEGs <- get(load("../inputfiles/LUAD_dataDEGs_deconvolution.rda")) 


#countshow many genese defined as both a tumoor supprosr and an onco gene by moonlight

Dualcount <- length(intersect(DataTS$x,DataOC$x))
TScount <- length(DataTS$x)
OCcount <- length(DataOC$x)
totalgene <- nrow(dataDEGs)






GeneDriverBarPlot=barplot(c(TScount,OCcount,totalgene),
        main = "Driver genes predicted by moonlight ",
        ylab = "Gene count",
        xlab = "Gene classification",
        ylim = c(0,14000),
        names.arg = c("Tumor suppressors","Onco genes","Genes entry"),
        col = c("green","red","gray"),
        horiz = FALSE,
        cex.axis = 1.5, 
        cex.lab=1.5)


dev.off()


