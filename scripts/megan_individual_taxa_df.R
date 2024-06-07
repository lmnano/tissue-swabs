megan_individual_taxa = function(inList){
  
  library(dplyr)
  library(tidyr)
  
  inListDF = bind_rows(inList) %>%
    pivot_wider(names_from = taxon, values_from = count)
  
  inListDF = as.data.frame(apply(inListDF, c(1,2), gsub, pattern = ",", replacement = ""))
  inListDF = as.data.frame(apply(inListDF, c(1,2), gsub, pattern = " ", replacement = ""))
  
  inListDF2 = apply(inListDF[,2:ncol(inListDF)], 2, as.integer)

  inListDFsum = rowSums(inListDF2, na.rm = TRUE)
  
  inListDF3 = apply(inListDF2, 2, function(x, sums){
    out = round(x / sums, 2)
    
  }, sums = inListDFsum)
  
  inListDFout = bind_cols(inListDF$name, inListDF3)
  
  colsNames = gsub(pattern = "other entries", replacement = "Other", x = colnames(inListDF3))
  
  colnames(inListDFout) = c("Sample", colsNames)
  
  inListDFout2 = inListDFout %>%
    mutate(Sample_type = grepl(pattern = "_s$", x = Sample)) %>%
    mutate(Sample_type = if_else(Sample_type == TRUE, "swab", "tissue"))
  
  return(inListDFout2)
  
}