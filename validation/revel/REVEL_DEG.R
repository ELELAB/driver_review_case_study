#Join of REVEL scores

#load tidyverse
library(tidyverse)

#open files, DEG does not have a header
REVEL_table=read_csv("./revel_with_transcript_ids")
DEG_anotatoin=read_csv(file="../input/cscape_input_converted_MAF_DEG.csv",col_names = FALSE)

#join DEG to REVEL scores based on chromosome, postion, ref and mutation
score_DEG=left_join(DEG_anotatoin,REVEL_table,by=c("X1"="chr","X2"="hg19_pos","X3"="ref","X4"="alt"))

#save files
write.csv(score_DEG,file="./scores_DEG")



