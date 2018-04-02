delvoom <- function(gDict, rawCountTable, deCoeff, design, molType, quantSoft, normMat){
  if(quantSoft == "kallisto"){
    normMat <- normMat/exp(rowMeans(log(normMat)))
    o <- log(calcNormFactors((rawCountTable)/normMat)) + log(colSums((rawCountTable)/normMat))
    y <- DGEList(rawCountTable)
    y$offset <- t(t(log(normMat)) + o)
  } else {
    
    y <- DGEList(counts = rawCountTable)
    y <- calcNormFactors(y)
  }
  
  v <- voom(y, design, plot = F)
  vfit <- lmFit(v, design)
  efit <- eBayes(vfit)
  
  toSyl <- topTable(efit, coef = deCoeff, n = Inf)
  
  if (molType == "gene"){
    toSyl <- toSyl[rownames(toSyl) %in% gDict$V2, ]
    rownames(toSyl) <- gDict[match(rownames(toSyl), gDict$V2), 1]
  }
  
  bySignif <- toSyl[order(log2(toSyl$adj.P.Val) * sign(toSyl$logFC)), ]
  return(bySignif)
}
