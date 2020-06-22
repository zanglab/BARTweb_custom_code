#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab 
#SBATCH -o out_convertID_refseq_append.log


# time for i in deg_FC1.5/*.txt;do convertIDs.pl $i human refseq > deg_FC1.5_refseq/$(basename $i);done

mkdir deg_FC1.5_refseq_append

time for i in deg_FC1.5_append/*.txt;do convertIDs.pl $i human refseq > deg_FC1.5_refseq_append/$(basename $i);done
