library(TFEA.ChIP)

indir="figure5_knockTF_validation/f1_deg_identification/deg_FC1.5_Entrez/"
outdir="tfea_deg_FC1.5"
dir.create(outdir)

# read all input gene list
allfiles<- list.files(path=indir,pattern=glob2rx("*txt"))
for (i in 1:length(allfiles)){
  genes <-read.table(file.path(indir,allfiles[i]),stringsAsFactors=FALSE)
  cont_genes <-contingency_matrix(genes[,1])
  results <- getCMstats( cont_genes ) 
  write.csv(results,file = file.path(outdir,allfiles[i]),row.names=TRUE)
}
