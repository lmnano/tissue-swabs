
tissue_swab_blast_top_hits_datatable = function(data, top_hits = 1){
  require(data.table)
  
  df = fread(data)
  
  blastCols = c("query_id", "subject_id", "perc_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "subject_sci_names", "subject_tax_ids", "subject_title")
  
  names(df) = blastCols
  
  d = data.table(df, key="query_id")
  d = d[, head(.SD, top_hits), by=query_id]
  
  return(d)
  
}
