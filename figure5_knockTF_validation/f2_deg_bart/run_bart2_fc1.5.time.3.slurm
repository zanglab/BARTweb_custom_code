#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_bart2_fc1.5.time.3.log

 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_178_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_178_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_02_73_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_02_73_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF2_DataSet_01_78_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF2_DataSet_01_78_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_117_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_117_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_118_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_118_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_133_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_133_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_134_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_134_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF1_DataSet_01_138_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF1_DataSet_01_138_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF2_DataSet_01_137_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF2_DataSet_01_137_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF2_DataSet_01_100_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF2_DataSet_01_100_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF4_DataSet_01_262_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF4_DataSet_01_262_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUND_DataSet_02_16_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUND_DataSet_02_16_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_164_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_164_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_172_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_172_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_230_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_291_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_291_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KAT2B_DataSet_02_184_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KAT2B_DataSet_02_184_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF3_DataSet_01_221_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF3_DataSet_01_221_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_339_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_339_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_341_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_341_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_340_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_340_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_342_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_342_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF9_DataSet_01_222_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF9_DataSet_01_222_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LDB1_DataSet_02_35_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LDB1_DataSet_02_35_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNA_DataSet_02_175_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNA_DataSet_02_175_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNB1_DataSet_02_195_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNB1_DataSet_02_195_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFG_DataSet_02_189_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFG_DataSet_02_189_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFK_DataSet_02_183_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFK_DataSet_02_183_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAF_DataSet_01_165_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAF_DataSet_01_165_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAX_DataSet_02_192_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAX_DataSet_02_192_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAZ_DataSet_02_199_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAZ_DataSet_02_199_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEF2A_DataSet_01_265_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEF2A_DataSet_01_265_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEIS2_DataSet_02_30_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEIS2_DataSet_02_30_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_207_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_207_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_24_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_24_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_25_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_25_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_79_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_79_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_02_85_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_02_85_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MSX1_DataSet_01_281_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MSX1_DataSet_01_281_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXD3_DataSet_02_83_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXD3_DataSet_02_83_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXI1_DataSet_02_173_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXI1_DataSet_02_173_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYBL2_DataSet_01_148_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_102_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_102_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_44_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_44_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_48_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_48_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_160_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_160_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_296_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_296_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_231_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_231_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_237_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_237_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_324_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_324_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NANOG_DataSet_01_94_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC1_DataSet_02_22_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC1_DataSet_02_22_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC3_DataSet_01_166_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC3_DataSet_01_166_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L1_DataSet_02_32_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L1_DataSet_02_32_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_335_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_335_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_60_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_60_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_02_174_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_02_174_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFKB1_DataSet_01_167_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFKB1_DataSet_01_167_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_131_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_131_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_287_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_288_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_288_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_355_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_355_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_02_76_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_02_76_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_01_238_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_01_238_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_02_24_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_02_24_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NKX2-1_DataSet_01_67_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NKX2-1_DataSet_01_67_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_149_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_149_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_65_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_65_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_93_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_93_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2C2_DataSet_02_77_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2C2_DataSet_02_77_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_07_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_07_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_89_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_89_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_90_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_90_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_02_19_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_02_19_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR3C1_DataSet_01_141_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR3C1_DataSet_01_141_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_01_174_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_01_174_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_02_84_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_02_84_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR5A2_DataSet_01_182_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR5A2_DataSet_01_182_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NRF1_DataSet_02_62_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NRF1_DataSet_02_62_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/OVOL2_DataSet_01_313_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/OVOL2_DataSet_01_313_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX2_DataSet_01_72_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX2_DataSet_01_72_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_26_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_27_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_27_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_359_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_359_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POLR3A_DataSet_01_139_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POLR3A_DataSet_01_139_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_05_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_05_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_42_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_42_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_95_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_95_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRDM14_DataSet_01_45_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRDM14_DataSet_01_45_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PROX1_DataSet_01_06_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PROX1_DataSet_01_06_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRPF4_DataSet_02_128_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRPF4_DataSet_02_128_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_158_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_158_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_99_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_99_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RAD21_DataSet_02_190_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RAD21_DataSet_02_190_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RARA_DataSet_01_51_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RARA_DataSet_01_51_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RCOR1_DataSet_02_204_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RCOR1_DataSet_02_204_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_30_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_30_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_80_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_80_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_02_39_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_02_39_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELB_DataSet_01_294_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_145_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_211_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_211_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_337_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_337_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_338_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
