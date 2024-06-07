
library(data.table)
source("tissue_swab_top_hits.R")

outFile = tissue_swab_blast_top_hits_datatable(data = "/path/to/in/file.tsv")
fwrite(outFile, "/path/to/output/outFile.csv")
