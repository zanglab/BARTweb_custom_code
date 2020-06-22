library(jsonlite)
									
indir="chea3_deg_FC1.5"
# read all input gene list
allfiles<- list.files(path=indir,pattern=glob2rx("*json"))
for (i in 1:length(allfiles)){
  results <- read_json(file.path(indir,allfiles[i]), simplifyVector = TRUE)
  print(allfiles[i])
  results_csv <- data.frame(
    TF = results$"Integrated--meanRank"$"TF",
    Score = results$"Integrated--meanRank"$"Score",
    Pvalue = results$"Integrated--topRank"$"Score"
  )
  write.csv(results_csv,file = strsplit(file.path(indir,allfiles[i]),".json")[[1]],row.names=TRUE)
}




