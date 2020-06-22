#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o slurm_deg_fc1.5/HMBOX1_DataSet_02_58_DnGenes.out 

findMotifs.pl ../f1_deg_identification/deg_FC1.5/HMBOX1_DataSet_02_58_DnGenes.txt human homer_deg_FC1.5/HMBOX1_DataSet_02_58_DnGenes -p 4
