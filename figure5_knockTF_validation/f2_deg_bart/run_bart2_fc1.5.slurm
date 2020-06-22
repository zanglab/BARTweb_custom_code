#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_bart2_fc1.5.log


bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_162_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_162_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_95_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO1_DataSet_02_95_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO2_DataSet_02_132_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AGO2_DataSet_02_132_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AHR_DataSet_01_283_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AHR_DataSet_01_283_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_108_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_108_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_109_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_109_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_345_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1A_DataSet_01_345_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1B_DataSet_01_298_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID1B_DataSet_01_298_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_193_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_193_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_50_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARID3A_DataSet_02_50_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_240_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_240_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_315_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_315_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ARNT_DataSet_01_317_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_192_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_278_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_278_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_323_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/AR_DataSet_01_323_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF2_DataSet_01_206_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF2_DataSet_01_206_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_01_33_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_01_33_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_180_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_180_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_54_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ATF3_DataSet_02_54_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_01_58_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_01_58_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_01_59_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_17_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_17_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_188_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BACH1_DataSet_02_188_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCL11A_DataSet_01_09_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCL11A_DataSet_01_09_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_01_214_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_112_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_112_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_142_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BCLAF1_DataSet_02_142_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BHLHE40_DataSet_02_186_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BHLHE40_DataSet_02_186_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRCA1_DataSet_02_202_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRCA1_DataSet_02_202_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_119_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_119_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_311_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/BRD4_DataSet_01_311_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CBFB_DataSet_01_43_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CBFB_DataSet_01_43_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CDX2_DataSet_01_357_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CDX2_DataSet_01_357_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPD_DataSet_01_74_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPD_DataSet_01_74_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_01_75_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_01_75_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_05_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_05_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_118_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_118_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_176_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CEBPZ_DataSet_02_176_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD1_DataSet_02_187_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD1_DataSet_02_187_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD2_DataSet_02_196_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CHD2_DataSet_02_196_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_04_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_04_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_194_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_194_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_195_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CREB1_DataSet_01_195_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTBP1_DataSet_02_41_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTBP1_DataSet_02_41_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTCF_DataSet_02_185_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/CTCF_DataSet_02_185_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DLX1_DataSet_02_55_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DLX1_DataSet_02_55_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DPF2_DataSet_02_14_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DPF2_DataSet_02_14_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_01_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_01_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_145_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/DROSHA_DataSet_02_145_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F4_DataSet_02_200_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F4_DataSet_02_200_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F6_DataSet_02_205_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/E2F6_DataSet_02_205_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_242_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_242_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_243_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR2_DataSet_01_243_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_212_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_212_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_34_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EGR3_DataSet_01_34_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EHF_DataSet_01_147_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EHF_DataSet_01_147_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_162_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_162_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_354_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK1_DataSet_01_354_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK3_DataSet_01_318_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ELK3_DataSet_01_318_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERF_DataSet_02_34_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERF_DataSet_02_34_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_17_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_17_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_216_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_216_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_28_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_28_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_316_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ERG_DataSet_01_316_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESR1_DataSet_01_01_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESR1_DataSet_01_01_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_115_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_115_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_280_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ESRRA_DataSet_01_280_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETS1_DataSet_01_244_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETS1_DataSet_01_76_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETS1_DataSet_01_76_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_161_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_161_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_179_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_179_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_36_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_36_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV4_DataSet_01_356_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV4_DataSet_01_356_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV5_DataSet_01_304_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV5_DataSet_01_304_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_122_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_122_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_167_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_167_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_12_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_12_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_13_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_13_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_14_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_14_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_15_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_15_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_16_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_16_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOSL1_DataSet_01_171_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOSL1_DataSet_01_171_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOS_DataSet_01_99_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOS_DataSet_01_99_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_104_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_104_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_105_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_233_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_233_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXK2_DataSet_02_31_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXK2_DataSet_02_31_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_122_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_122_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_223_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_223_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_77_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_77_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_02_37_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_02_37_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_106_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_106_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_270_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_270_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_321_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_321_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_322_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_322_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO3_DataSet_01_267_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO3_DataSet_01_267_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_208_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_208_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_209_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_209_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_210_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_210_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_234_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_234_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_88_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_88_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_124_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_124_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_155_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_155_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_113_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_113_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_190_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_191_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA1_DataSet_02_60_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA1_DataSet_02_60_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_193_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_193_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_61_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_61_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_02_28_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_02_28_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_120_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_120_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_163_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_163_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_347_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_347_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA6_DataSet_01_247_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA6_DataSet_01_247_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GLI2_DataSet_01_143_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GRHL2_DataSet_01_289_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GRHL3_DataSet_01_336_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GRHL3_DataSet_01_336_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2B_DataSet_01_362_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2B_DataSet_01_362_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_110_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_110_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_157_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_157_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDAC8_DataSet_02_51_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDAC8_DataSet_02_51_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDGF_DataSet_02_104_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDGF_DataSet_02_104_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HES1_DataSet_01_266_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HES1_DataSet_01_266_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIC1_DataSet_01_205_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIC1_DataSet_01_205_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_224_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_224_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_348_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_348_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_86_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_86_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_87_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_87_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HINFP_DataSet_02_36_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HINFP_DataSet_02_36_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMBOX1_DataSet_02_58_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMBOX1_DataSet_02_58_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGA1_DataSet_02_75_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGA1_DataSet_02_75_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN1_DataSet_02_79_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN1_DataSet_02_79_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN3_DataSet_02_45_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN3_DataSet_02_45_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1A_DataSet_01_02_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1A_DataSet_01_02_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1B_DataSet_01_103_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1B_DataSet_01_103_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPLL_DataSet_02_125_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPLL_DataSet_02_125_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_160_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_160_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_96_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_96_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXA10_DataSet_01_358_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXA1_DataSet_01_110_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXA1_DataSet_01_110_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB4_DataSet_02_64_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB4_DataSet_02_64_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_169_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_169_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_170_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_170_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXC6_DataSet_01_333_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXC6_DataSet_01_333_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXD9_DataSet_01_70_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXD9_DataSet_01_70_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_112_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_112_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_155_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_155_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_178_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_178_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_02_73_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_02_73_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF2_DataSet_01_78_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF2_DataSet_01_78_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_117_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_117_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_118_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_118_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_133_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_133_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_134_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ID1_DataSet_01_134_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF1_DataSet_01_138_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF1_DataSet_01_138_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF2_DataSet_01_137_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IKZF2_DataSet_01_137_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF2_DataSet_01_100_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF2_DataSet_01_100_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF4_DataSet_01_262_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/IRF4_DataSet_01_262_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUND_DataSet_02_16_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUND_DataSet_02_16_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_164_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_164_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_172_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_172_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_230_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_291_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/JUN_DataSet_01_291_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KAT2B_DataSet_02_184_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KAT2B_DataSet_02_184_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF3_DataSet_01_221_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF3_DataSet_01_221_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_339_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_339_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_341_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF4_DataSet_01_341_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_340_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_340_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_342_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF7_DataSet_01_342_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF9_DataSet_01_222_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/KLF9_DataSet_01_222_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LDB1_DataSet_02_35_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LDB1_DataSet_02_35_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNA_DataSet_02_175_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNA_DataSet_02_175_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNB1_DataSet_02_195_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/LMNB1_DataSet_02_195_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFG_DataSet_02_189_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFG_DataSet_02_189_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFK_DataSet_02_183_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAFK_DataSet_02_183_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAF_DataSet_01_165_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAF_DataSet_01_165_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAX_DataSet_02_192_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAX_DataSet_02_192_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAZ_DataSet_02_199_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MAZ_DataSet_02_199_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEF2A_DataSet_01_265_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEF2A_DataSet_01_265_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEIS2_DataSet_02_30_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MEIS2_DataSet_02_30_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_207_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_207_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_24_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_24_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_25_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_25_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_79_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_01_79_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_02_85_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MITF_DataSet_02_85_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MSX1_DataSet_01_281_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MSX1_DataSet_01_281_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXD3_DataSet_02_83_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXD3_DataSet_02_83_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXI1_DataSet_02_173_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MXI1_DataSet_02_173_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYBL2_DataSet_01_148_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_102_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_102_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_44_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_44_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_48_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYB_DataSet_01_48_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_160_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_160_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_296_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYCN_DataSet_01_296_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_231_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_231_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_237_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_237_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_324_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/MYC_DataSet_01_324_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NANOG_DataSet_01_94_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC1_DataSet_02_22_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC1_DataSet_02_22_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC3_DataSet_01_166_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFATC3_DataSet_01_166_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L1_DataSet_02_32_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L1_DataSet_02_32_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_335_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_335_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_60_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_01_60_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_02_174_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFE2L2_DataSet_02_174_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFKB1_DataSet_01_167_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFKB1_DataSet_01_167_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_131_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_131_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_287_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_288_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_288_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_355_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_01_355_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_02_76_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYA_DataSet_02_76_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_01_238_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_01_238_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_02_24_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NFYB_DataSet_02_24_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NKX2-1_DataSet_01_67_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NKX2-1_DataSet_01_67_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_149_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_149_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_65_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_65_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_93_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NONO_DataSet_02_93_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2C2_DataSet_02_77_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2C2_DataSet_02_77_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_07_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_07_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_89_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_89_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_90_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_01_90_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_02_19_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR2F2_DataSet_02_19_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR3C1_DataSet_01_141_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR3C1_DataSet_01_141_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_01_174_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_01_174_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_02_84_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR4A1_DataSet_02_84_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR5A2_DataSet_01_182_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NR5A2_DataSet_01_182_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NRF1_DataSet_02_62_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/NRF1_DataSet_02_62_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/OVOL2_DataSet_01_313_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/OVOL2_DataSet_01_313_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX2_DataSet_01_72_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX2_DataSet_01_72_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_26_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_27_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_27_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_359_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PAX3_DataSet_01_359_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POLR3A_DataSet_01_139_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POLR3A_DataSet_01_139_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_05_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_05_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_42_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_42_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_95_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_95_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRDM14_DataSet_01_45_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRDM14_DataSet_01_45_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PROX1_DataSet_01_06_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PROX1_DataSet_01_06_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRPF4_DataSet_02_128_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PRPF4_DataSet_02_128_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_158_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_158_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_99_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/PSIP1_DataSet_02_99_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RAD21_DataSet_02_190_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RAD21_DataSet_02_190_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RARA_DataSet_01_51_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RARA_DataSet_01_51_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RCOR1_DataSet_02_204_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RCOR1_DataSet_02_204_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_30_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_30_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_80_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_01_80_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_02_39_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELA_DataSet_02_39_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RELB_DataSet_01_294_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_145_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_211_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_211_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_337_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_337_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_338_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/REST_DataSet_01_338_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RFX5_DataSet_02_61_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RFX5_DataSet_02_61_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RNF2_DataSet_02_68_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RNF2_DataSet_02_68_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_159_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_159_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_22_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_22_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_314_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_314_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_47_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_47_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_93_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX1_DataSet_01_93_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX2_DataSet_01_62_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/RUNX2_DataSet_01_62_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_277_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_277_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_98_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SALL4_DataSet_01_98_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SETDB1_DataSet_02_52_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SETDB1_DataSet_02_52_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_146_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_146_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_271_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SF1_DataSet_01_271_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX4_DataSet_01_269_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX4_DataSet_01_269_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX5_DataSet_02_78_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SIX5_DataSet_02_78_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD2_DataSet_01_351_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD3_DataSet_01_249_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD3_DataSet_01_352_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD5_DataSet_02_48_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMAD5_DataSet_02_48_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMARCA4_DataSet_02_47_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SMARCA4_DataSet_02_47_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI1_DataSet_01_69_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI1_DataSet_01_69_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_225_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_225_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_292_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNAI2_DataSet_01_292_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNRNP70_DataSet_02_126_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SNRNP70_DataSet_02_126_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_300_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_300_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_360_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX17_DataSet_01_360_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_349_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_349_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_350_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_97_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX2_DataSet_01_97_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_140_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_140_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_260_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX4_DataSet_01_260_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_308_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_308_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_309_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SOX9_DataSet_01_309_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_01_81_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_01_81_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_02_67_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP1_DataSet_02_67_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP2_DataSet_02_42_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP2_DataSet_02_42_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP3_DataSet_01_290_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SP3_DataSet_01_290_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_107_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_107_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_135_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SPDEF_DataSet_01_135_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRF_DataSet_02_21_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRF_DataSet_02_21_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRSF3_DataSet_02_127_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/SRSF3_DataSet_02_127_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_154_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_154_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_21_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_21_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_83_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_01_83_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_02_23_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT1_DataSet_02_23_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT2_DataSet_02_18_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT2_DataSet_02_18_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_168_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_168_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_180_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_184_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_184_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_334_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_334_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_49_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_49_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_84_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT3_DataSet_01_84_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_182_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_182_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_72_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5A_DataSet_02_72_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_239_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_239_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_279_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT5_DataSet_01_279_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_50_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_50_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_85_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_01_85_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_02_57_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/STAT6_DataSet_02_57_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_106_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_106_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_140_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_140_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_163_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAF15_DataSet_02_163_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_152_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_152_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_217_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_217_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_40_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_41_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_01_41_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_02_201_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TAL1_DataSet_02_201_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_116_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_116_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_129_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_129_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_169_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_169_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_80_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TARDBP_DataSet_02_80_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_178_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_178_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_43_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TBL1XR1_DataSet_02_43_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF4_DataSet_01_185_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF4_DataSet_01_185_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF7L1_DataSet_01_346_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TCF7L1_DataSet_01_346_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD2_DataSet_02_33_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD2_DataSet_02_33_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_01_151_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_01_151_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_02_40_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TEAD4_DataSet_02_40_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_150_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_150_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_312_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_312_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_53_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP2C_DataSet_01_53_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFAP4_DataSet_01_295_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFDP1_DataSet_02_71_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TFDP1_DataSet_02_71_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/THAP1_DataSet_02_46_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/THAP1_DataSet_02_46_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TLX1_DataSet_01_246_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TLX1_DataSet_01_263_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_10_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_10_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_11_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_11_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_175_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_175_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_286_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_286_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_31_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_31_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_37_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_37_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_38_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_38_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_68_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP53_DataSet_01_68_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_183_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_183_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_196_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_196_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_197_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_197_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_198_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_198_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_199_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_199_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_200_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_200_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_201_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_201_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_202_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_202_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_325_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_325_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_327_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_327_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_328_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_328_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_330_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_330_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_39_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TP63_DataSet_01_39_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TRIM28_DataSet_02_63_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/TRIM28_DataSet_02_63_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_165_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_165_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_97_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/U2AF1_DataSet_02_97_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_198_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_198_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_56_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/UBTF_DataSet_02_56_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF1_DataSet_02_26_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF1_DataSet_02_26_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF2_DataSet_02_25_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/USF2_DataSet_02_25_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/VEZF1_DataSet_02_179_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/VEZF1_DataSet_02_179_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WHSC1_DataSet_02_66_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WHSC1_DataSet_02_66_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WT1_DataSet_01_08_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/WT1_DataSet_01_08_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_203_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_203_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_204_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/XBP1_DataSet_01_204_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_255_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_255_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_257_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_257_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_258_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_258_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX1_DataSet_01_259_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_153_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_153_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_98_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YBX3_DataSet_02_98_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YY1_DataSet_01_18_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/YY1_DataSet_01_18_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZBTB33_DataSet_02_82_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZBTB33_DataSet_02_82_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZC3H11A_DataSet_02_44_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZC3H11A_DataSet_02_44_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB1_DataSet_01_232_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB2_DataSet_01_56_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB2_DataSet_01_56_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZEB2_DataSet_01_57_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_01_144_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_01_144_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_01_157_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_172_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_172_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_206_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZFX_DataSet_02_206_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZIC2_DataSet_01_121_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZIC2_DataSet_01_121_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZMIZ1_DataSet_02_181_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZMIZ1_DataSet_02_181_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF100_DataSet_01_124_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF100_DataSet_01_124_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF143_DataSet_02_59_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF143_DataSet_02_59_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF180_DataSet_02_177_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF180_DataSet_02_177_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF24_DataSet_01_353_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF24_DataSet_01_353_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF254_DataSet_01_127_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF254_DataSet_01_127_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF25_DataSet_01_319_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF25_DataSet_01_319_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_203_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_203_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_27_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF384_DataSet_02_27_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF395_DataSet_02_69_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF395_DataSet_02_69_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF431_DataSet_01_128_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF431_DataSet_01_128_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF708_DataSet_01_130_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF708_DataSet_01_130_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF85_DataSet_01_126_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5
wait

bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ZNF85_DataSet_01_126_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5 &
