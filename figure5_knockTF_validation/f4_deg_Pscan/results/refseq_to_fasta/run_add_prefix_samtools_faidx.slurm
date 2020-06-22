#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab

#Run program

fai="Homo_sapiens_1000up_0down.fasta.fai.uniq"
fasta_file="figure5_knockTF_validation/f4_deg_Pscan/results/Human_and_Mouse/Homo_sapiens_1000up_0down.fasta"

# fa file from DEG using 1.5FC
indir="figure5_knockTF_validation/f1_deg_identification/deg_FC1.5_refseq"
outdir="deg_FC1.5_refseq_fa"
mkdir ${outdir}
for i in ${indir}/*.txt
do 
    cat $i|awk '{print "hg38_refGene_"$0}' > ${outdir}/$(basename $i)
    cat $fai ${outdir}/$(basename $i) |sort |uniq -d > ${outdir}/$(basename $i).fai 
    xargs samtools faidx ${fasta_file} < ${outdir}/$(basename $i).fai > ${outdir}/$(basename $i).fai.fa
done


