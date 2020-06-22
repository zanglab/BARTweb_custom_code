#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_bart2_fc1.5.time.1.log


time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_162_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_162_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_95_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_95_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO2_DataSet_02_132_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO2_DataSet_02_132_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AHR_DataSet_01_283_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AHR_DataSet_01_283_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_108_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_108_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_109_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_109_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_345_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_345_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1B_DataSet_01_298_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1B_DataSet_01_298_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_193_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_193_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_50_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_50_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_240_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_240_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_315_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_315_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_317_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_192_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_278_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_278_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_323_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_323_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF2_DataSet_01_206_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF2_DataSet_01_206_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_01_33_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_01_33_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_180_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_180_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_54_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_54_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_01_58_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_01_58_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_01_59_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_17_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_17_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_188_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_188_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCL11A_DataSet_01_09_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCL11A_DataSet_01_09_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_01_214_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_112_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_112_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_142_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_142_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BHLHE40_DataSet_02_186_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BHLHE40_DataSet_02_186_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRCA1_DataSet_02_202_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRCA1_DataSet_02_202_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_119_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_119_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_311_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_311_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CBFB_DataSet_01_43_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CBFB_DataSet_01_43_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CDX2_DataSet_01_357_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CDX2_DataSet_01_357_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPD_DataSet_01_74_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPD_DataSet_01_74_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_01_75_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_01_75_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_05_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_05_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_118_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_118_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_176_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_176_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD1_DataSet_02_187_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD1_DataSet_02_187_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD2_DataSet_02_196_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD2_DataSet_02_196_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_04_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_04_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_194_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_194_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_195_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_195_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTBP1_DataSet_02_41_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTBP1_DataSet_02_41_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTCF_DataSet_02_185_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTCF_DataSet_02_185_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DLX1_DataSet_02_55_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DLX1_DataSet_02_55_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DPF2_DataSet_02_14_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DPF2_DataSet_02_14_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_01_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_01_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_145_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_145_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F4_DataSet_02_200_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F4_DataSet_02_200_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F6_DataSet_02_205_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F6_DataSet_02_205_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_242_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_242_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_243_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_243_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_212_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_212_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_34_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_34_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EHF_DataSet_01_147_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EHF_DataSet_01_147_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_162_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_162_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_354_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_354_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK3_DataSet_01_318_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK3_DataSet_01_318_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERF_DataSet_02_34_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERF_DataSet_02_34_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_17_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_17_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_216_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_216_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_28_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_28_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_316_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_316_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESR1_DataSet_01_01_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESR1_DataSet_01_01_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_115_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_115_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_280_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_280_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETS1_DataSet_01_244_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETS1_DataSet_01_76_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETS1_DataSet_01_76_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 
