deEdgeR <- function(gDict, rawCountTable, deCoeff, design, molType, quantSoft, normMat){
  if(quantSoft == "kallisto"){
    normMat <- normMat/exp(rowMeans(log(normMat)))
    o <- log(calcNormFactors((rawCountTable)/normMat)) + log(colSums((rawCountTable)/normMat))
    y <- DGEList(rawCountTable)
    y$offset <- t(t(log(normMat)) + o)
    sy <- estimateDisp(y, design)
    fit <- glmFit(y, design, sy$tagwise.dispersion)
    
  } else {
    
    y <- DGEList(counts = rawCountTable)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
  }
  
  toSyl <- topTags(glmLRT(fit, coef = deCoeff), n = Inf)$table
  if (molType == "gene"){
    toSyl <- toSyl[rownames(toSyl) %in% gDict$V2, ]
    rownames(toSyl) <- gDict[match(rownames(toSyl), gDict$V2), 1]
  }
  
  bySignif <- toSyl[order(log2(toSyl$FDR) * sign(toSyl$logFC)), ]
  return(bySignif)
}
