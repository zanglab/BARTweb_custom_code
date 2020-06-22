#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o out_convertID_Entrez_append.log


# mkdir deg_FC1.5_Entrez
# 
# time for i in deg_FC1.5/*.txt;do convertIDs.pl $i human gene > deg_FC1.5_Entrez/$(basename $i);done

mkdir deg_FC1.5_Entrez_append

time for i in deg_FC1.5_append/*.txt;do convertIDs.pl $i human gene > deg_FC1.5_Entrez_append/$(basename $i);done

