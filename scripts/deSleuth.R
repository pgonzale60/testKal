deSleuth <- function(gDict, kal_dirs, slTest, molType, titlNote, sampleInfo, toModMat){
  s2c <- dplyr::mutate(sampleInfo, sample = Run_s, path = kal_dirs)
  
  t2g <- dplyr::rename(gDict, target_id = V1, ens_gene = V2)
  
  if (molType == "gene"){
    so <- sleuth_prep(s2c, target_mapping = t2g,
                      aggregation_column = 'ens_gene',
                      extra_bootstrap_summary = TRUE)
    
  } else {
    so <- sleuth_prep(s2c, target_mapping = t2g,
                      extra_bootstrap_summary = TRUE)
  }
  
  
  
  
  so <- sleuth_fit(so, eval(parse(text=toModMat)), 'full')
  so <- sleuth_wt(so, which_beta = slTest)
  cont.res <- sleuth_results(so, slTest, show_all = F)
  
  if (molType == "gene"){
    bySignif <- cont.res[order(log2(cont.res$qval) * sign(cont.res$b)), ] %>%
      dplyr::mutate(target_id = gDict[match(target_id, gDict$V2), 1]) %>%
      dplyr::select(target_id, b, qval)
    
  } else {
    bySignif <- cont.res[order(log2(cont.res$qval) * sign(cont.res$b)), ] %>%
      dplyr::select(target_id, b, qval)
  }
  
  outFile <- paste("../output/DE/", molType, "/sleuth/kallisto_",
                   srp, "_", titlNote, ".tsv.gz", sep = "")
  write.table(bySignif, file = gzfile(outFile),
              row.names=F, col.names=F, sep="\t", quote=FALSE)
}
