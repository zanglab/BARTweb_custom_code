library(httr)
library(jsonlite)
									
indir="figure5_knockTF_validation/f1_deg_identification/deg_FC1.5/"
outdir="chea3_deg_FC1.5"
dir.create(outdir)

# read all input gene list
allfiles<- list.files(path=indir,pattern=glob2rx("*txt"))
for (i in 1:length(allfiles)){

genes <-read.table(file.path(indir,allfiles[i]),stringsAsFactors=FALSE)
url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes[,1])
response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")
results = fromJSON(json)
write_json(results,file.path(outdir,paste(allfiles[i],".json",sep="")))

results_csv <- data.frame(
#  Rank = results$"Integrated--meanRank"$"Rank",
  TF = results$"Integrated--meanRank"$"TF",
  Score = results$"Integrated--meanRank"$"Score"
#   Overlapping_Genes = results$"Integrated--meanRank"$"Overlapping_Genes",
)

write.csv(results_csv,file = file.path(outdir,allfiles[i]),row.names=TRUE)

}
