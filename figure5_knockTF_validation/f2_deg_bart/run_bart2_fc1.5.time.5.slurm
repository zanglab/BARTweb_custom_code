#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_bart2_fc1.5.time.log

 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD2_DataSet_02_33_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD2_DataSet_02_33_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_01_151_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_01_151_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_02_40_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_02_40_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_150_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_150_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_312_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_312_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_53_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_53_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP4_DataSet_01_295_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFDP1_DataSet_02_71_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFDP1_DataSet_02_71_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/THAP1_DataSet_02_46_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/THAP1_DataSet_02_46_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TLX1_DataSet_01_246_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TLX1_DataSet_01_263_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_10_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_10_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_11_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_11_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_175_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_175_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_286_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_286_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_31_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_31_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_37_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_37_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_38_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_38_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_68_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_68_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_183_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_183_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_196_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_196_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_197_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_197_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_198_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_198_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_199_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_199_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_200_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_200_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_201_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_201_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_202_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_202_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_325_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_325_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_327_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_327_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_328_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_328_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_330_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_330_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_39_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_39_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TRIM28_DataSet_02_63_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TRIM28_DataSet_02_63_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_165_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_165_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_97_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_97_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_198_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_198_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_56_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_56_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF1_DataSet_02_26_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF1_DataSet_02_26_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF2_DataSet_02_25_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF2_DataSet_02_25_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/VEZF1_DataSet_02_179_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/VEZF1_DataSet_02_179_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WHSC1_DataSet_02_66_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WHSC1_DataSet_02_66_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WT1_DataSet_01_08_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WT1_DataSet_01_08_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_203_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_203_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_204_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_204_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_255_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_255_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_257_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_257_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_258_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_258_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_259_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_153_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_153_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_98_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_98_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YY1_DataSet_01_18_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YY1_DataSet_01_18_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZBTB33_DataSet_02_82_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZBTB33_DataSet_02_82_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZC3H11A_DataSet_02_44_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZC3H11A_DataSet_02_44_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB1_DataSet_01_232_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB2_DataSet_01_56_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB2_DataSet_01_56_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB2_DataSet_01_57_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_01_144_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_01_144_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_01_157_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_172_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_172_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_206_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_206_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZIC2_DataSet_01_121_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZIC2_DataSet_01_121_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZMIZ1_DataSet_02_181_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZMIZ1_DataSet_02_181_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF100_DataSet_01_124_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF100_DataSet_01_124_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF143_DataSet_02_59_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF143_DataSet_02_59_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF180_DataSet_02_177_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF180_DataSet_02_177_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF24_DataSet_01_353_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF24_DataSet_01_353_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF254_DataSet_01_127_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF254_DataSet_01_127_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF25_DataSet_01_319_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF25_DataSet_01_319_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_203_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_203_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_27_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_27_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF395_DataSet_02_69_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF395_DataSet_02_69_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF431_DataSet_01_128_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF431_DataSet_01_128_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF708_DataSet_01_130_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF708_DataSet_01_130_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF85_DataSet_01_126_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF85_DataSet_01_126_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
