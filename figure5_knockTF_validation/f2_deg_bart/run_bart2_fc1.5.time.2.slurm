#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=100000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_bart2_fc1.5.time.2.log



time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_161_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_161_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_179_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_179_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_36_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV1_DataSet_01_36_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV4_DataSet_01_356_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV4_DataSet_01_356_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV5_DataSet_01_304_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/ETV5_DataSet_01_304_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_122_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_122_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_167_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/EWSR1_DataSet_02_167_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_12_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_12_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_13_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_13_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_14_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_14_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_15_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_15_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_16_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FLI1_DataSet_01_16_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOSL1_DataSet_01_171_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOSL1_DataSet_01_171_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOS_DataSet_01_99_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOS_DataSet_01_99_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_104_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_104_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_105_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_233_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXA1_DataSet_01_233_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXK2_DataSet_02_31_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXK2_DataSet_02_31_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_122_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_122_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_223_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_223_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_77_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_01_77_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_02_37_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXM1_DataSet_02_37_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_106_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_106_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_270_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_270_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_321_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_321_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_322_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO1_DataSet_01_322_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO3_DataSet_01_267_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXO3_DataSet_01_267_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_208_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_208_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_209_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_209_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_210_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_210_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_234_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_234_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_88_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FOXP1_DataSet_01_88_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_124_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_124_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_155_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/FUS_DataSet_02_155_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_113_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_113_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_190_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GABPA_DataSet_01_191_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA1_DataSet_02_60_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA1_DataSet_02_60_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_193_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_193_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_61_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_01_61_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_02_28_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA2_DataSet_02_28_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_120_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_120_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_163_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_163_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_347_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA3_DataSet_01_347_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA6_DataSet_01_247_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GATA6_DataSet_01_247_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GLI2_DataSet_01_143_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GRHL2_DataSet_01_289_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GRHL3_DataSet_01_336_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GRHL3_DataSet_01_336_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2B_DataSet_01_362_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2B_DataSet_01_362_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_110_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_110_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_157_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/GTF2F1_DataSet_02_157_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDAC8_DataSet_02_51_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDAC8_DataSet_02_51_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDGF_DataSet_02_104_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HDGF_DataSet_02_104_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HES1_DataSet_01_266_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HES1_DataSet_01_266_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIC1_DataSet_01_205_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIC1_DataSet_01_205_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_224_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_224_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_348_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_348_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_86_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_86_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_87_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HIF1A_DataSet_01_87_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HINFP_DataSet_02_36_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HINFP_DataSet_02_36_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMBOX1_DataSet_02_58_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMBOX1_DataSet_02_58_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGA1_DataSet_02_75_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGA1_DataSet_02_75_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN1_DataSet_02_79_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN1_DataSet_02_79_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN3_DataSet_02_45_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HMGN3_DataSet_02_45_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1A_DataSet_01_02_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1A_DataSet_01_02_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1B_DataSet_01_103_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNF1B_DataSet_01_103_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPLL_DataSet_02_125_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPLL_DataSet_02_125_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_160_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_160_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_96_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HNRNPL_DataSet_02_96_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXA10_DataSet_01_358_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXA1_DataSet_01_110_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXA1_DataSet_01_110_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB4_DataSet_02_64_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB4_DataSet_02_64_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_169_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_169_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_170_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXB7_DataSet_01_170_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXC6_DataSet_01_333_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXC6_DataSet_01_333_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXD9_DataSet_01_70_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
 

time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HOXD9_DataSet_01_70_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_112_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_112_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_155_DnGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time 
time bart2 geneset -i ../f1_deg_identification/deg_FC1.5/HSF1_DataSet_01_155_UpGenes.txt -s hg38 --outdir bart2_deg_FC1.5.time
