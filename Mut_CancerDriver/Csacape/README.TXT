run and analysys of CScape for MAF inner joined by input genes for Moonlight.


Here Cscape:somatic's table of driver mutations prediciton is used on MAF file, which have been fillted to only contain genes used in prior analyse for Moonlight. 

two MAF file was process by:
1. inner joined by gene found in DEG and a second with oncogenic mediatos from Moonlight.
2. lifted from CGhr 38 to CGHr 38.
3. Colums extraction on chromosome posion, mutaionpostion, ref and mutated.

see conversion_MAF_DEG.R and conversion_MAF_DEG.R
 

CScape's table was the runed by the bash scrips.


Analyse of Cscape resutles were perforemd in analyse_CScape.R
