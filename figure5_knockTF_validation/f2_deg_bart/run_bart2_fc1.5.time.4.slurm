#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_bart2_fc1.5.time.4.log

 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_338_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RFX5_DataSet_02_61_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RFX5_DataSet_02_61_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RNF2_DataSet_02_68_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RNF2_DataSet_02_68_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_159_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_159_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_22_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_22_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_314_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_314_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_47_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_47_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_93_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_93_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX2_DataSet_01_62_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX2_DataSet_01_62_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_277_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_277_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_98_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_98_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SETDB1_DataSet_02_52_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SETDB1_DataSet_02_52_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_146_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_146_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_271_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_271_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX4_DataSet_01_269_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX4_DataSet_01_269_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX5_DataSet_02_78_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX5_DataSet_02_78_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD2_DataSet_01_351_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD3_DataSet_01_249_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD3_DataSet_01_352_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD5_DataSet_02_48_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD5_DataSet_02_48_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMARCA4_DataSet_02_47_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMARCA4_DataSet_02_47_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI1_DataSet_01_69_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI1_DataSet_01_69_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_225_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_225_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_292_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_292_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNRNP70_DataSet_02_126_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNRNP70_DataSet_02_126_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_300_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_300_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_360_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_360_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_349_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_349_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_350_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_97_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_97_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_140_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_140_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_260_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_260_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_308_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_308_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_309_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_309_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_01_81_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_01_81_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_02_67_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_02_67_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP2_DataSet_02_42_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP2_DataSet_02_42_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP3_DataSet_01_290_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP3_DataSet_01_290_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_107_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_107_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_135_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_135_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRF_DataSet_02_21_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRF_DataSet_02_21_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRSF3_DataSet_02_127_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRSF3_DataSet_02_127_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_154_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_154_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_21_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_21_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_83_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_83_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_02_23_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_02_23_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT2_DataSet_02_18_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT2_DataSet_02_18_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_168_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_168_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_180_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_184_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_184_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_334_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_334_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_49_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_49_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_84_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_84_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_182_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_182_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_72_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_72_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_239_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_239_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_279_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_279_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_50_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_50_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_85_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_85_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_02_57_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_02_57_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_106_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_106_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_140_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_140_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_163_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_163_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_152_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_152_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_217_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_217_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_40_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_41_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_41_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_02_201_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_02_201_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_116_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_116_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_129_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_129_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_169_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_169_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_80_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_80_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_178_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_178_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_43_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_43_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF4_DataSet_01_185_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF4_DataSet_01_185_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF7L1_DataSet_01_346_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF7L1_DataSet_01_346_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 
