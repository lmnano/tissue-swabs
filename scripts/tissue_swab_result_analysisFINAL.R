library(Rsamtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)

### Comparison of RAD Loci generated from minimally invasive swab samples and invasive tissue samples ###
## RAD loci generating in STACKS

### Comparison of different mapping algorithms ###
# Mapped data not uploaded due to size and available upon request

###TABLE S2###
#Path to data
bamPath = "./data/tissue_swabs_comparisons/bowtie"
blastPath = "./data/tissue_swabs_comparisons/blast"
#Getting names for bam and index files and importing bam files along with their indexes into a list
bamFiles = list.files(bamPath, pattern = ".*bam", full.names = TRUE)
bamIndex = list.files(bamPath, pattern = ".*bai", full.names = TRUE)
bamIndex = gsub(pattern = "()\\.bai", replacement = "\\1", x = bamIndex)
bamNames = list.files(bamPath, pattern = ".*bam")
bamNames = gsub(pattern = "()\\.bam", replacement = "\\1", x = bamNames)
bamImport = lapply(1:length(bamFiles), function(x, b, i){
  out = scanBam(b[x], index = i[x])
},b = bamFiles, i = bamIndex)
names(bamImport) = bamNames
#Extracting query and subject sequence names from bam list
querySubjectBowtie = lapply(1:length(bamImport), function(x, bam_list, comboNames){
  qname = bam_list[[x]][[1]]$qname
  rname = as.character(bam_list[[x]][[1]]$rname)
  df = data.frame(qname, rname)
  df$qname = as.numeric(df$qname)
  isNumber = sum(!is.na(as.numeric(df$rname)))
  if (isNumber > 0) {
    df$rname = as.numeric(df$rname)
  }else{
    pos = as.character(bam_list[[x]][[1]]$pos)
    df$rname = paste(df$rname, pos, sep = "_")
  }
  df$rname = na_if(as.character(df$rname), "NA_NA")
  df$unique_subject = !df$rname %in% df$rname[duplicated(df$rname)]
  colnames(df) = c("query_id", paste("subject", comboNames[[x]], sep = "_"), paste("unique_subject", comboNames[[x]], sep = "_"))
  df
}, bam_list = bamImport, comboNames = bamNames)
names(querySubjectBowtie) = bamNames
#Reading blast results and taking query and subject columns
blastFiles = list.files(blastPath, pattern = "blast.*csv", full.names = TRUE)
blastNames = list.files(blastPath, pattern = "blast.*csv")
blastNames = gsub(pattern = "()_Top.*", replacement = "\\1", x = blastNames)
blastTopResults = lapply(blastFiles, read.csv, header = TRUE)
names(blastTopResults) = blastNames
querySubjectBlast = lapply(blastTopResults, function(x){
  if (class(x$subject_id) == "character") {
    x$subject_id = paste(x$subject_id, as.character(x$s_start), sep ="_")
  }
  out = x[,1:2]
  out$unique_subject_blast = !out$subject_id %in% out$subject_id[duplicated(out$subject_id)]
  out
})
#Joining results into tables by different combinations: whole tissue vs. swabs, whole swabs vs tissue, split swabs and split tissue vs. whole swabs and tissue and split swabs and tissue vs split swabs and tissue. Six combinations in total. Genome data added in the same manner, four extra tables, 10 in total.
querySubjectAllSvT = querySubjectBowtie$bowtie2AllSvT_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2AllSvT_local_sorted) %>%
  left_join(querySubjectBlast$blastAllSvT)
querySubjectAllSvT$unique_subject_bowtie2AllSvT_global_sorted[is.na(querySubjectAllSvT$subject_bowtie2AllSvT_global_sorted)] = NA
querySubjectAllSvT$unique_subject_bowtie2AllSvT_local_sorted[is.na(querySubjectAllSvT$subject_bowtie2AllSvT_local_sorted)] = NA
querySubjectAllTvS = querySubjectBowtie$bowtie2AllTvS_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2AllTvS_local_sorted) %>%
  left_join(querySubjectBlast$blastAllTvS)
querySubjectAllTvS$unique_subject_bowtie2AllTvS_global_sorted[is.na(querySubjectAllTvS$subject_bowtie2AllTvS_global_sorted)] = NA
querySubjectAllTvS$unique_subject_bowtie2AllTvS_local_sorted[is.na(querySubjectAllTvS$subject_bowtie2AllTvS_local_sorted)] = NA
querySubjectSplitAllSvT = querySubjectBowtie$bowtie2SplitAllSvT_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2SplitAllSvT_local_sorted) %>%
  left_join(querySubjectBlast$blastSplitAllSvT)
querySubjectSplitAllSvT$unique_subject_bowtie2SplitAllSvT_global_sorted[is.na(querySubjectSplitAllSvT$subject_bowtie2SplitAllSvT_global_sorted)] = NA
querySubjectSplitAllSvT$unique_subject_bowtie2SplitAllSvT_local_sorted[is.na(querySubjectSplitAllSvT$subject_bowtie2SplitAllSvT_local_sorted)] = NA
querySubjectSplitAllTvS = querySubjectBowtie$bowtie2SplitAllTvS_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2SplitAllTvS_local_sorted) %>%
  left_join(querySubjectBlast$blastSplitAllTvS)
querySubjectSplitAllTvS$unique_subject_bowtie2SplitAllTvS_global_sorted[is.na(querySubjectSplitAllTvS$subject_bowtie2SplitAllTvS_global_sorted)] = NA
querySubjectSplitAllTvS$unique_subject_bowtie2SplitAllTvS_local_sorted[is.na(querySubjectSplitAllTvS$subject_bowtie2SplitAllTvS_local_sorted)] = NA
querySubjectSplitSvT = querySubjectBowtie$bowtie2SplitSvT_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2SplitSvT_local_sorted) %>%
  left_join(querySubjectBlast$blastSplitSvT)
querySubjectSplitSvT$unique_subject_bowtie2SplitSvT_global_sorted[is.na(querySubjectSplitSvT$subject_bowtie2SplitSvT_global_sorted)] = NA
querySubjectSplitSvT$unique_subject_bowtie2SplitSvT_local_sorted[is.na(querySubjectSplitSvT$subject_bowtie2SplitSvT_local_sorted)] = NA
querySubjectSplitTvS = querySubjectBowtie$bowtie2SplitTvS_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2SplitTvS_local_sorted) %>%
  left_join(querySubjectBlast$blastSplitTvS)
querySubjectSplitTvS$unique_subject_bowtie2SplitTvS_global_sorted[is.na(querySubjectSplitTvS$subject_bowtie2SplitTvS_global_sorted)] = NA
querySubjectSplitTvS$unique_subject_bowtie2SplitTvS_local_sorted[is.na(querySubjectSplitTvS$subject_bowtie2SplitTvS_local_sorted)] = NA
querySubjectAllSvG = querySubjectBowtie$bowtie2AllSvG_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2AllSvG_local_sorted) %>%
  left_join(querySubjectBlast$blastAllSvG)
querySubjectAllSvG$unique_subject_bowtie2AllSvG_global_sorted[is.na(querySubjectAllSvG$subject_bowtie2AllSvG_global_sorted)] = NA
querySubjectAllSvG$unique_subject_bowtie2AllSvG_local_sorted[is.na(querySubjectAllSvG$subject_bowtie2AllSvG_local_sorted)] = NA
querySubjectAllTvG = querySubjectBowtie$bowtie2AllTvG_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2AllTvG_local_sorted) %>%
  left_join(querySubjectBlast$blastAllTvG)
querySubjectAllTvG$unique_subject_bowtie2AllTvG_global_sorted[is.na(querySubjectAllTvG$subject_bowtie2AllTvG_global_sorted)] = NA
querySubjectAllTvG$unique_subject_bowtie2AllTvG_local_sorted[is.na(querySubjectAllTvG$subject_bowtie2AllTvG_local_sorted)] = NA
querySubjectSplitSvG = querySubjectBowtie$bowtie2SplitSvG_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2SplitSvG_local_sorted) %>%
  left_join(querySubjectBlast$blastSplitSvG)
querySubjectSplitSvG$unique_subject_bowtie2SplitSvG_global_sorted[is.na(querySubjectSplitSvG$subject_bowtie2SplitSvG_global_sorted)] = NA
querySubjectSplitSvG$unique_subject_bowtie2SplitSvG_local_sorted[is.na(querySubjectSplitSvG$subject_bowtie2SplitSvG_local_sorted)] = NA
querySubjectSplitTvG = querySubjectBowtie$bowtie2SplitTvG_global_sorted %>%
  left_join(querySubjectBowtie$bowtie2SplitTvG_local_sorted) %>%
  left_join(querySubjectBlast$blastSplitTvG)
querySubjectSplitTvG$unique_subject_bowtie2SplitTvG_global_sorted[is.na(querySubjectSplitTvG$subject_bowtie2SplitTvG_global_sorted)] = NA
querySubjectSplitTvG$unique_subject_bowtie2SplitTvG_local_sorted[is.na(querySubjectSplitTvG$subject_bowtie2SplitTvG_local_sorted)] = NA
#Basic metrics and graphs
querySubjectList = list(querySubjectAllSvT, querySubjectAllTvS, querySubjectSplitAllSvT, querySubjectSplitAllTvS, querySubjectSplitSvT, querySubjectSplitTvS, querySubjectAllSvG, querySubjectAllTvG, querySubjectSplitSvG, querySubjectSplitTvG)
names(querySubjectList) = c("querySubjectAllSvT", "querySubjectAllTvS", "querySubjectSplitAllSvT", "querySubjectSplitAllTvS", "querySubjectSplitSvT", "querySubjectSplitTvS", "querySubjectAllSvG", "querySubjectAllTvG", "querySubjectSplitSvG", "querySubjectSplitTvG")
querySubjectList = lapply(querySubjectList, function(x){
  x$bt2_global_local_same = x[,2] == x[,4] | (is.na(x[,2]) & is.na(x[,4]))
  x$bt2_global_local_same[is.na(x$bt2_global_local_same)] = FALSE
  x$bt2_global_blast_same = x[,2] == x[,6] | (is.na(x[,2]) & is.na(x[,6]))
  x$bt2_global_blast_same[is.na(x$bt2_global_blast_same)] = FALSE
  x$bt2_local_blast_same = x[,4] == x[,6] | (is.na(x[,4]) & is.na(x[,6]))
  x$bt2_local_blast_same[is.na(x$bt2_local_blast_same)] = FALSE
  x$all_same = x[,2] == x[,4] & x[,2] == x[,6] & x[,4] == x[,6] | (is.na(x[,2]) & is.na(x[,4]) & is.na(x[,6]))
  x$all_same[is.na(x$all_same)] = FALSE
  x
})
querySubjectList_basicMetrics = lapply(querySubjectList, function(x){
  out = data.frame(Bt2_perc_global_uniq_withNA = sum(x[,3], na.rm = TRUE)/nrow(x))
  out$Bt2_perc_local_uniq_withNA = sum(x[,5], na.rm = TRUE)/nrow(x)
  out$blast_perc_uniq_withNA = sum(x[,7], na.rm = TRUE)/nrow(x)
  out$Bt2_perc_global_uniq_woNA = sum(x[,3], na.rm = TRUE)/(nrow(x)-sum(is.na(x[,3])))
  out$Bt2_perc_local_uniq_woNA = sum(x[,5], na.rm = TRUE)/(nrow(x)-sum(is.na(x[,5])))
  out$blast_perc_uniq_woNA = sum(x[,7], na.rm = TRUE)/(nrow(x)-sum(is.na(x[,7])))
  out$Bt2_perc_hits_global = 1-(sum(is.na(x[,2]))/nrow(x))
  out$Bt2_perc_hits_global = 1-(sum(is.na(x[,2]))/nrow(x))
  out$Bt2_perc_hits_local = 1-(sum(is.na(x[,4]))/nrow(x))
  out$blast_perc_hits = 1-(sum(is.na(x[,6]))/nrow(x))
  out$perc_Bt2_global_local_same = sum(x[,8])/nrow(x)
  out$perc_Bt2_global_blast_same = sum(x[,9])/nrow(x)
  out$perc_Bt2_local_blast_same = sum(x[,10])/nrow(x)
  out$perc_all_same = sum(x[,11])/nrow(x)
  out$Bt2_perc_local_NOTuniq_withNA = sum(x[,5] == FALSE, na.rm = TRUE)/nrow(x)
  out$Bt2_perc_local_NA_withNA = sum(is.na(x[,5]))/nrow(x)
  out$Bt2_perc_global_uniq_NUM = sum(x[,3], na.rm = TRUE)
  out$Bt2_perc_global_NUM = nrow(x) - sum(is.na(x[,3]))
  out$Bt2_perc_local_uniq_NUM = sum(x[,5], na.rm = TRUE)
  out$Bt2_perc_local_NUM = nrow(x) - sum(is.na(x[,5]))
  out$blast_perc_uniq_NUM = sum(x[,7], na.rm = TRUE)
  out$blast_perc_NUM = nrow(x) - sum(is.na(x[,7]))
  out
})
querySubjectMetrics = do.call(rbind, querySubjectList_basicMetrics)
querySubjectMetricsWITH_NUMS = querySubjectMetrics
querySubjectMetrics = querySubjectMetrics[,c(1:15)]
rownames(querySubjectMetrics) = gsub("querySubject(.*)", replacement = "\\1", x = rownames(querySubjectMetrics))
querySubjectMetrics$names = rownames(querySubjectMetrics)
querySubjectMetricsFINAL = querySubjectMetrics %>%
  filter(names %in% c("AllSvT", "SplitAllSvT", "SplitSvT", "AllSvG")) %>%
  mutate(names = gsub(pattern = "^AllSvT$", replacement = "swab vs. tissue full", x = names)) %>%
  mutate(names = gsub(pattern = "^SplitAllSvT$", replacement = "swab split vs. tissue", x = names)) %>%
  mutate(names = gsub(pattern = "^SplitSvT$", replacement = "swab split vs. tissue split", x = names)) %>%
  mutate(names = gsub(pattern = "^AllSvG$", replacement = "swabs vs. genome", x = names)) %>%
  t()
querySubjectMetricsFINAL = as.data.frame(querySubjectMetricsFINAL)
colnames(querySubjectMetricsFINAL) = querySubjectMetricsFINAL["names",]
querySubjectMetricsFINAL = querySubjectMetricsFINAL[1:(nrow(querySubjectMetricsFINAL)-3),]
querySubjectMetricsFINAL = apply(querySubjectMetricsFINAL, 2, as.numeric)
querySubjectMetricsFINAL = apply(querySubjectMetricsFINAL, 2, round, 3)
rownames(querySubjectMetricsFINAL) = c("unique with NAs bowtie gl.", "unique with NAs bowtie lo.", "unique with NAs blast", "unique without NAs bowtie gl.", "unique without NAs bowtie lo.", "unique without NAs blast", "total bowtie gl.", "total bowtie lo.", "total blast", "overlapping bowtie gl. vs lo.", "overlapping bowtie gl. vs blast", "overlapping lo. vs blast", "all similar")

###FIGURE S2A###
table_perc_hits = pivot_longer(querySubjectMetrics[, c(1:13, 16)], cols = c(7:9), names_to = "metric", values_to = "values") %>%
  select(11:13)
table_perc_hitsFINAL = table_perc_hits %>%
  filter(names %in% c("AllSvT", "SplitSvT", "SplitAllSvT")) %>%
  mutate(metric = (gsub(pattern = "^Bt2_perc_hits_global$", replacement = "Bowtie2 global", x = metric))) %>%
  mutate(metric = (gsub(pattern = "^Bt2_perc_hits_local$", replacement = "Bowtie2 local", x = metric))) %>%
  mutate(metric = (gsub(pattern = "^blast_perc_hits$", replacement = "Blast", x = metric)))
ggplot(table_perc_hitsFINAL, aes(x = names, y = values, fill = metric)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylim(0:1) +
  labs(title = "Percentage of found sequences", x = "Sequence concatenation", y = "Proportion of mapped loci")

###FIGURE 4###
querySubjectMetrics_bestMethodTG = querySubjectMetrics %>%
  dplyr::filter(names == "AllSvT" | names == "AllSvG") %>%
  pivot_longer(cols = c(1:9, 14,15), names_to = "metric", values_to = "values") %>%
  select(c("names", "metric", "values")) %>%
  dplyr::filter(grepl(pattern = "local", x = metric)) %>%
  mutate(prettyMetric = metric) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_uniq_withNA", replacement = "unique matches", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_uniq_woNA", replacement = "Unique percentage without unmtched", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_hits_local", replacement = "Percentage All", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_NOTuniq_withNA", replacement = "non-unique matches", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_NA_withNA", replacement = "unmatched", x = prettyMetric)) %>%
  group_by(names) %>%
  arrange(values, .by_group = TRUE) %>%
  mutate(prettyMetric = fct_inorder(factor(prettyMetric)))
querySubjectMetrics_bestMethodSTG = querySubjectMetrics %>%
  dplyr::filter( names == "AllSvG" | names == "AllTvG") %>%
  pivot_longer(cols = c(1:9, 14,15), names_to = "metric", values_to = "values") %>%
  select(c("names", "metric", "values")) %>%
  dplyr::filter(grepl(pattern = "local", x = metric)) %>%
  mutate(prettyMetric = metric) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_uniq_withNA", replacement = "unique matches", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_uniq_woNA", replacement = "Unique percentage without unmtched", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_hits_local", replacement = "Percentage All", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_NOTuniq_withNA", replacement = "non-unique matches", x = prettyMetric)) %>%
  mutate(prettyMetric = gsub(pattern = "Bt2_perc_local_NA_withNA", replacement = "unmatched", x = prettyMetric)) %>%
  group_by(names) %>%
  arrange(values, .by_group = TRUE) %>%
  mutate(prettyMetric = fct_inorder(factor(prettyMetric)))
querySubjectMetrics_bestMethodFINAL = querySubjectMetrics_bestMethodSTG %>%
  filter(names == "AllTvG") %>%
  bind_rows(querySubjectMetrics_bestMethodTG) %>%
  mutate(names = (gsub(pattern = "^AllTvG$", replacement = "tissue catalogue vs. genome", x = names))) %>%
  mutate(names = (gsub(pattern = "^AllSvG$", replacement = "swab catalogue vs. genome", x = names))) %>%
  mutate(names = (gsub(pattern = "^AllSvT$", replacement = "swab vs. tissue catalogue", x = names)))
querySubjectMetrics_bestMethodFINAL %>%
  filter(grepl(".*withNA$", metric)) %>%
  ggplot(aes(x = names, y = values, fill = prettyMetric)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #scale_y_continuous(n.breaks=6) +
  ylim(0:1) +
  labs(y = "Proportion")

###Exogenous DNA results anlysis###

###FIGURE 5###
taxaTablesFullNames = list.files(path = "./data/megan_output_limit/", full.names = TRUE)
taxaTablesNames = list.files(path = "./data/megan_output_limit/")
taxaTableList = lapply(taxaTablesFullNames, read.table, sep = "\t")
names(taxaTableList) = taxaTablesNames
intxtdf = lapply(taxaTablesFullNames, function(x){
  out = readLines(x)
  out = out[6:(length(out) - 1)]
  outdf = as.data.frame(out)
  outdf = as.data.frame(str_split_fixed(outdf$out, pattern = ": ", 2))
  outdf$V1 = str_squish(outdf$V1)
  colnames(outdf) = c("taxon", "count")
  return(outdf)
})
names(intxtdf) = taxaTablesNames
intxtdf2 = lapply(1:length(intxtdf), function(x, indata){
  df = indata[[x]]
  namedf = names(indata[x])
  domain = grepl(pattern = "bacteria", x = namedf)
  if (domain == TRUE) {
    out = df %>%
      mutate(count = gsub(pattern = ",", replacement = "", x = count)) %>%
      mutate(count = as.integer(count)) %>%
      mutate(perc = count / sum(count)) %>%
      mutate(taxon2 = if_else(perc < 0.05, "Others", taxon)) %>%
      #mutate(taxon2 = if_else(taxon == "unclassified Bacteria", "Unclassified", taxon2)) %>%
      group_by(taxon2) %>%
      summarise(percSum = sum(perc), .groups = "keep") %>%
      ungroup() %>%
      arrange((percSum)) %>%
      mutate(fctlevel = 1:nrow(.)) %>%
      mutate(fctlevel = if_else(taxon2 == "Others", as.integer(nrow(.) + 10000), fctlevel)) %>%
      arrange(fctlevel) %>%
      mutate(taxon2 = factor(taxon2)) %>%
      mutate(taxon2 = fct_inorder(taxon2))
  }else{
    out = df %>%
      mutate(count = gsub(pattern = ",", replacement = "", x = count)) %>%
      mutate(count = as.integer(count)) %>%
      mutate(perc = count / sum(count)) %>%
      mutate(taxon2 = taxon) %>%
      #mutate(taxon2 = if_else(taxon == "unclassified Bacteria", "Unclassified", taxon2)) %>%
      group_by(taxon2) %>%
      summarise(percSum = sum(perc), .groups = "keep") %>%
      ungroup() %>%
      arrange((percSum)) %>%
      mutate(fctlevel = 1:nrow(.)) %>%
      mutate(fctlevel = if_else(taxon2 == "Others", as.integer(nrow(.) + 10000), fctlevel)) %>%
      arrange(fctlevel) %>%
      mutate(taxon2 = factor(taxon2)) %>%
      mutate(taxon2 = fct_inorder(taxon2))
  }
  return(out)
}, indata = intxtdf)
names(intxtdf2) = taxaTablesName
lapply(1:length(intxtdf2), function(x, indata){
  df = indata[[x]]
  dfname = names(indata[x])
  taxon = gsub(pattern = ".*_outfmt0_(.*)_summary.txt", replacement = "\\1", x = dfname)
  ggplot(df, aes(x="", y=percSum, fill = fct_rev(taxon2))) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    labs(title = taxon,
         x = "taxon",
         y = "Percentage") +
    guides(fill = guide_legend(title = taxon, reverse = TRUE)) +
    theme_void()
  #ggsave(filename = paste0("megan_", dfname, taxon, "_newBLASTrun_limit", ".eps"), path = "./figures")
}, indata = intxtdf2)

###TABLE S5###
topSpecies = function(taxonTable, topn = 10){
  library(dplyr)
  topSpecies = taxonTable %>%
    ungroup() %>%
    group_by(species) %>%
    summarise(number = n()) %>%
    arrange(desc(number)) %>%
    slice_max(number, n = topn)
  return(topSpecies)
}
taxaTableUnfiltered = taxaTableList[[4]]
taxaTableUnfilteredName = taxaTablesNames[[4]]
comboName = sub(pattern  = "taxa_.*(blast.*)_.*", replacement = "\\1", x = taxaTableUnfilteredName)
taxaTable = taxaTableUnfiltered %>%
  group_by(query_id) %>%
  arrange(evalue, desc(perc_identity), .by_group = TRUE) %>%
  slice_head(n = 1)
topSpeciesAll_Swab_vs_genome = topSpecies(taxaTable)
taxaTableUnfiltered = taxaTableList[[6]]
taxaTableUnfilteredName = taxaTablesNames[[6]]
comboName = sub(pattern  = "taxa_.*(blast.*)_.*", replacement = "\\1", x = taxaTableUnfilteredName)
taxaTable = taxaTableUnfiltered %>%
  group_by(query_id) %>%
  arrange(evalue, desc(perc_identity), .by_group = TRUE) %>%
  slice_head(n = 1)
topSpeciesAll_Swab_vs_tissue = topSpecies(taxaTable)
taxaTableUnfiltered = taxaTableList[[8]]
taxaTableUnfilteredName = taxaTablesNames[[8]]
comboName = sub(pattern  = "taxa_.*(blast.*)_.*", replacement = "\\1", x = taxaTableUnfilteredName)
taxaTable = taxaTableUnfiltered %>%
  group_by(query_id) %>%
  arrange(evalue, desc(perc_identity), .by_group = TRUE) %>%
  slice_head(n = 1)
topSpeciesAll_Tissue_vs_genome = topSpecies(taxaTable)

