sigDeStat <- function(srp, colToRelevl, ref, toModMat, titlNote, molType,
                      quantSoft, statSoft, deCoeff, slTest,
                      dictF = "../data/transcriptomes/pcMus_trans_gene.tsv",
                      sleuthF){
  library("edgeR")
  library("tximport")
  library("sleuth")
  
  
  if(!(quantSoft %in% c("featureCounts", "kallisto", "stringtie"))){
    message("Quantification software only allows: featureCounts, kallisto or stringtie")
    stop()
  }
  
  if(!(statSoft %in% c("edgeR", "sleuth", "limma-voom"))){
    message("Quantification software only allows: edgeR, sleuth or limma-voom")
    stop()
  }
  
  if(!(molType %in% c("gene", "trans", "utrans"))){
    message("Quantification software only allows: gene, trans or utrans")
    stop()
  }
  
  if(quantSoft == "featureCounts" & molType != "gene"){
    message("featureCounts only works for gene counts")
    stop()
  }
  
  if(quantSoft == "stringtie" & molType == "utrans"){
    message("stringtie did not estimated counts per representative transcript (utrans).")
    stop()
  }
  
  if(statSoft == "sleuth" & quantSoft != "kallisto"){
    message("Sleuth only is compatible with kallisto results.")
    stop()
  }
  
  if(!file.exists(sleuthF) & statSoft != "sleuth"){
    message("Sleuth result not found. Need this to filter any other strategy.")
    stop()
  }
  
  
  # Read sample info
  sampleInfo <- read.table(paste("../data/metadata/", srp, "_metadat.tsv", sep = ""), header=TRUE, sep="\t", as.is = T)
  sampleInfo[, colToRelevl] <- relevel(factor(sampleInfo[, colToRelevl]), ref=ref)
  
  # Read dictionary
  gDict <- read.table(dictF, header=F, sep="\t", as.is = T)
  
  
  # Sleuth rutine
  if(quantSoft == "kallisto" & statSoft == "sleuth"){
    if(molType == "utrans")
      kalT <- "uKallisto"
    else
      kalT <- "cKallisto"
    
    kal_dirs <- list.dirs(paste("../output/counts/", kalT, "/", srp, "/", sep = ""),
                          recursive = F, full.names = T)
    
    deSleuth(gDict, kal_dirs, slTest, molType, titlNote, sampleInfo, toModMat)
    
  } else {
    # Read sleuth accepted features
    gKal <- read.table(sleuthF, header=F, sep="\t", as.is = T)
    # And set dictionary according to molecule type
    if(molType == "gene"){
      kalGens <- unique(gDict[match(gKal$V1, gDict$V1), 2])
    } else{
      kalGens <- gKal$V1
    }
    
    # Create design matrix
    design <- model.matrix(eval(parse(text=toModMat)), data=sampleInfo)
    
    # Initialize empty variable which might be filled by tximport
    normMat <- 0
    
    #### Read counts ####
    
    # From featureCounts
    if(quantSoft == "featureCounts"){
      rawCountTable <- read.table(paste("../output/counts/featureCounts/",
                                        srp, "/", srp, ".counts.tsv",
                                        sep = ""), header=TRUE, sep="\t",
                                  row.names=1)[, -c(1:5), drop = F]
      
      colnames(rawCountTable) <- sub("(.+).sorted.bam", "\\1", colnames(rawCountTable))
    }
    
    # From stringTie
    if(quantSoft == "stringtie"){
      if(molType == "gene")
        stFile <- "gene_count_matrix.csv"
      else
        stFile <- "transcript_count_matrix.csv"
      
      rawCountTable <- read.table(paste("../output/counts/stringtie/",
                                        srp, "/", stFile,
                                        sep = ""), header=TRUE, sep=",",
                                  row.names=1)
    }
    
    # From kallisto
    if(quantSoft == "kallisto" & statSoft != "sleuth"){
      if(molType != "utrans")
        kalT <- "cKallisto"
      else
        kalT <- "uKallisto"
      
      files <- file.path(paste("../output/counts/", kalT, "/", srp, "/", sep = ""), sampleInfo$Run_s, "abundance.h5")
      names(files) <- sampleInfo$Run_s
      
      if(molType == "gene"){
        txOut <- FALSE
      }
      else{
        txOut <- TRUE
      }
      
      
      txi.kallisto <- tximport(files, type = "kallisto", txOut = txOut, tx2gene = gDict)
      rawCountTable <- txi.kallisto$counts
      normMat <- txi.kallisto$length[rownames(rawCountTable) %in% kalGens, ]
    }
    
    # Filter only accepted features by Kallisto
    kalAccp <- rownames(rawCountTable) %in% kalGens
    rawCountTable <- rawCountTable[kalAccp, ]
    
    #### Statistical analyses ####
    # EdgeR
    if(statSoft == "edgeR"){
      bySignif <- deEdgeR(gDict, rawCountTable, deCoeff, design, molType, quantSoft, normMat)
    }
    
    # Limma-voom
    if(statSoft == "limma-voom"){
      bySignif <- delvoom(gDict, rawCountTable, deCoeff, design, molType, quantSoft, normMat)
    }
    
    # Write results
    outFile <- paste("../output/DE/", molType, "/", statSoft, "/", quantSoft, "_",
                     srp, "_", titlNote, ".tsv.gz", sep = "")
    write.table(bySignif, file = gzfile(outFile),
                row.names=T, col.names=F, sep="\t", quote=FALSE)
  }
}
