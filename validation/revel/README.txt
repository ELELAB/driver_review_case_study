Here pathogenicity was predicted on REVEL for the LUAD data set.


REVEL's prediciton table was downlaode from https://sites.google.com/site/revelgenomics/downloads?authuser=0


The table was marged on the previees made MAF-like files which was used for CScape. 

The code for the margering is in REVEL_DEG.R and the bash script  for running on the server is run_REVEL_DEG
for marrgering tideverse inner.join fucntions was used.


The analysed was performed in analysis_REVEL.R
Since the oncogenic mediators predited is a subset on genes in DEG, the oncogenic mediatos set was found by inner.join oncogenic mediators genes
to the results of DEG, insted of running REVEL on the second maf files to avoid wait on run time. 

  
