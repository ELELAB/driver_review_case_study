###Moonlight Cancer  driver genes lung cancer subtype luad, subsubtype luad1

#Import Moonlight
library(MoonlightR)


#data load
dataDEGs <- get(load("../inputfiles/LUAD_dataDEGs_deconvolution.rda")) 
DataFilt <- get(load("../inputfiles/LUAD_dataFilt_adjusted.rda"))

#Find normal vs tumor tissue in samples
#"NT" for normal tissue, "TP" for primoary tumor
#library(TCGAbiolinks)
#barcodes <- colnames(DataFilt)
#TCGAquery_SampleTypes(barcodes, c("NB","NT","NBC","NEBV","NBM",))
#TCGAquery_SampleTypes(barcodes, c("TP"))
#clin <- GDCquery_clinic("TCGA-LUAD", type = "clinical", save.csv = TRUE)


#Enrichment( finds biological processes overrepresenteded in data)
dataFEA <- FEA(DEGsmatrix = dataDEGs)
#Save file
saveRDS(dataFEA,file = "../results/dataFEA.rds")


# Gene regulatory network used to idenfity genes closesly related to high deregulated gense
dataGRN <- GRN(TFs = rownames(dataDEGs), 
               DEGsmatrix = dataDEGs,
               DiffGenes = TRUE,
               normCounts = dataFilt)
saveRDS(dataGRN,file = "../results/dataGRN.rds")
#Save file


#selecet bioproees from the enrichemtn and
dataURA <- URA(dataGRN = dataGRN, 
               DEGsmatrix = dataDEGs, 
               BPname = c("apoptosis",
                          "proliferation of cells"),
               nCores=4)
#Save file
saveRDS(dataURA,file = "../results/dataURA.rds")


#findesDriver genese
dataPRA <- PRA(dataURA = dataURA, 
               BPname = c("apoptosis",
                          "proliferation of cells"),
               thres.role = 0)
#save file
saveRDS(dataPRA,file = "../results/dataPRA.rds")


#extra the TSG and OCG names to a list
CancerGenes <- list("TSG"=names(dataPRA$TSG), "OCG"=names(dataPRA$OCG))

#write the TSG and OCG into seperate files
write.csv(CancerGenes$TSG, file = "../results/Moon_CancerGenes_EB_TSG.csv", col.names = TRUE)  
write.csv(CancerGenes$OCG, file = "../results/Moon_CancerGenes_EB_OCG.csv", col.names = TRUE)  

