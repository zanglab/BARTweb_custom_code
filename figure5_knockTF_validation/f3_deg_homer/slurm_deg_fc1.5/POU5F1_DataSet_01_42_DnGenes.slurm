#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o slurm_deg_fc1.5/POU5F1_DataSet_01_42_DnGenes.out 

findMotifs.pl ../f1_deg_identification/deg_FC1.5/POU5F1_DataSet_01_42_DnGenes.txt human homer_deg_FC1.5/POU5F1_DataSet_01_42_DnGenes -p 4
