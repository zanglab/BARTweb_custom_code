#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o slurm_deg_fc1.5/AGO2_DataSet_02_132_UpGenes.out 

findMotifs.pl ../f1_deg_identification/deg_FC1.5/AGO2_DataSet_02_132_UpGenes.txt human homer_deg_FC1.5/AGO2_DataSet_02_132_UpGenes -p 4
