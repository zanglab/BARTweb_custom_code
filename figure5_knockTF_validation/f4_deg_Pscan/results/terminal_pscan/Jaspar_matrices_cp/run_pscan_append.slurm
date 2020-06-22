#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab


indir="figure5_knockTF_validation/f4_deg_Pscan/results/refseq_to_fasta/deg_FC1.5_refseq_fa_append"
for i in ${indir}/*.fa
do
    pscan -q ${i} -m Homo_sapiens_1000_0.short_matrix
done


