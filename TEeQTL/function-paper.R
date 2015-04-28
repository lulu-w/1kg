valCols <- function(x){
  m = x[, grep("id|type", colnames(x), invert = T, perl = T)]
  rownames(m) = x$id
  m
}

summaryPoly <- function(x){
    #Table 1
    m = valCols(x)
    rownames(m) = x$id
    res = x[which(rowSums(m)>0), "id"]
    res1 = gsub("_[0-9]+$", "", res, perl = T)
    table(res1)
}

getType <- function(x){
    df = data.frame(id = names(x), af = x)
    df$type = gsub("_[0-9]+$", "", names(x), perl = T)
    df$type = gsub("_umary_.*$", "", df$type, perl = T )
    df1 = df
    df1$type = "All"
    df = rbind(df, df1)
}

formatTable1 <- function(x){
  rownames(x) = gsub("_umary_.*$", "", rownames(x), perl = T )
  x = x[, c(group.labs, colnames(x)[grep("Union|All", colnames(x), perl = T)])]
}

getAF <- function(x){
    m = valCols(x)
    af = rowSums(m)/(ncol(m)*2)
    #af = rowSums(m) 
}

filterAF <- function(x, af.cut){
    m = valCols(x)
    rownames(m) = x$id
    af = rowSums(m)/(ncol(m)*2)
    x[which(af > af.cut), ]
    #af[which(af > af.cut)]
}


filterEQTL <- function(x, p.cut, fdr.cut){
    y = x[which(x$pvalue < p.cut & x$FDR < fdr.cut), ]
    y$snp.type = gsub("_[0-9]+$", "", y$snps, perl = T)
    y$snp.type = gsub("_umary_.*$", "", y$snp.type, perl = T )
    y
}


runME <- function( SNP_file_name, 
                   snps_location_file_name, 
                   expression_file_name, 
                   gene_location_file_name, 
                   covariates_file_name, 
                   useModel = modelLINEAR, 
                   par.plot = 100, 
                   verbose = FALSE, 
                   output_file_name = tempfile(), 
                   base.dir = ".", 
                   pvOutputThreshold = 1e-2, 
                   noFDRsaveMemory = FALSE)
{
  require(MatrixEQTL)
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  #useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  # Type of plots
  #pvalue.hist: #breaks(ie 100, 1000), "qqplot"
  
  errorCovariance = numeric();
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  ## Run the analysis
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = verbose,
    pvalue.hist = par.plot,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = noFDRsaveMemory);
  
  unlink(output_file_name);
  
  return(me)
}


runCisTrans <- function( SNP_file_name, 
                         snps_location_file_name, 
                         expression_file_name, 
                         gene_location_file_name, 
                         covariates_file_name, 
                         useModel = modelLINEAR, 
                         par.plot = 100, 
                         verbose = FALSE, 
                         output_file_name_cis = tempfile(), 
                         output_file_name_tra = tempfile(), 
                         base.dir = ".", 
                         pvOutputThreshold_cis = 2e-2, 
                         pvOutputThreshold_tra = 1e-2, 
                         cisDist = 1e6, 
                         noFDRsaveMemory = FALSE)
{
  require(MatrixEQTL)
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  #useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  # Type of plots
  #pvalue.hist: #breaks(ie 100, 1000), "qqplot"
  
  errorCovariance = numeric();  
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }

  
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
  
  
  ## Run the analysis

  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = verbose,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = par.plot,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  return(me)
}


myPlotQQ <- function (x, cex = 0.5, pch = 19, xlim = NULL, ylim = NULL, ...) 
{
  if (x$param$pvalue.hist == FALSE) {
    warning("Cannot plot p-value distribution: the information was not recorded.\nUse pvalue.hist!=FALSE.")
    return(invisible())
  }
  if (x$param$pvalue.hist == "qqplot") {
    xmin = 1/max(x$cis$ntests, x$all$ntests)
    ymax = NULL
    if (!is.null(ylim)) {
      ymax = ylim[2]
    }
    else {
      ymax = -log10(min(x$cis$eqtls$pvalue[1], x$cis$hist.bins[c(FALSE, 
                                                                 x$cis$hist.counts > 0)][1], x$all$eqtls$pvalue[1], 
                        x$all$hist.bins[c(FALSE, x$all$hist.counts > 
                                            0)][1], x$trans$eqtls$pvalue[1], x$trans$hist.bins[c(FALSE, 
                                                                                                 x$trans$hist.counts > 0)][1], na.rm = TRUE)) + 
        0.1
    }
    if (ymax == 0) {
      ymax = -log10(.Machine$double.xmin)
    }
    if (!is.null(ymax)) 
      ylim = c(0, ymax)
    if (is.null(xlim)) 
      xlim = c(0, -log10(xmin/1.5))
    plot(numeric(), numeric(), xlab = "-Log10(p-value), theoretical", 
         ylab = "-Log10(p-value), observed", xlim = xlim, 
         ylim = ylim, xaxs = "i", yaxs = "i", ...)
    lines(c(0, 1000), c(0, 1000), col = "gray")
    if ((x$param$pvOutputThreshold > 0) && (x$param$pvOutputThreshold.cis > 
                                              0)) {
      MatrixEQTL:::.qqme(x$cis, "red", cex, pch, ...)
      MatrixEQTL:::.qqme(x$trans, "blue", cex, pch, ...)
      title(paste("QQ-plot for", formatC(x$cis$ntests, 
                                         big.mark = ",", format = "f", digits = 0), "local and", 
                  formatC(x$trans$ntests, big.mark = ",", format = "f", 
                          digits = 0), "distant gene-SNP p-values"))
      lset = c(1, 2, 4)
    }
    else if (x$param$pvOutputThreshold.cis > 0) {
      MatrixEQTL:::.qqme(x$cis, "red", cex, pch, ...)
      title(paste("QQ-plot for", formatC(x$cis$ntests, 
                                         big.mark = ",", format = "f", digits = 0), "local gene-SNP p-values"))
      lset = c(1, 4)
    }
    else {
      MatrixEQTL:::.qqme(x$all, "blue", cex, pch, ...)
      title(paste("QQ-plot for all", formatC(x$all$ntests, 
                                             big.mark = ",", format = "f", digits = 0), "gene-SNP p-values"))
      lset = c(3, 4)
    }
    legend("top", c("Local p-values", "Distant p-values", 
                        "All p-values", "diagonal")[lset], col = c("red", 
                                                                   "blue", "blue", "gray")[lset], text.col = c("red", 
                                                                                                               "blue", "blue", "gray")[lset], pch = 20, lwd = 1, 
           pt.cex = c(1, 1, 1, 0)[lset])
  }
  else {
    if ((x$param$pvOutputThreshold > 0) && (x$param$pvOutputThreshold.cis > 
                                              0)) {
      par(mfrow = c(2, 1))
      .histme(x$cis, "", " local", ...)
      tran = list(hist.counts = x$all$hist.counts - x$cis$hist.counts, 
                  hist.bins = x$all$hist.bins, ntests = x$all$ntests - 
                    x$cis$ntests)
      .histme(x$trans, "", " distant", ...)
      par(mfrow = c(1, 1))
    }
    else if (x$param$pvOutputThreshold.cis > 0) {
      .histme(x$cis, "", " local", ...)
    }
    else {
      .histme(x$all, "all ", "", ...)
    }
  }
  return(invisible())
}

myChisq <- function(x, expected = covar2){
  y = unlist(x)
  x = ifelse(y > 0, 1, 0)
  chi.test.df = data.frame(id = names(x), presence = x, stringsAsFactors = F)
  mdf = merge(chi.test.df, expected, by.x = "id", by.y = "sample" )
  tb = table(mdf$group, mdf$presence)
  chisq.test(tb)$p.value
}


checkExp <- function(x, verbose = F){
  #lapply(1:5, function(x) ploteqtl(alleqtl[x, ]))
  #color points by "genotype" or "population"
  if(verbose) print(x)
  
  require(ggplot2)
  sid = as.character(x$snps)
  gid = as.character(x$gene)  
  pval = as.numeric(x$pvalue)
  
  xsnps <- snps[which(snps$id == sid ), ]
  xexp <- geneExp[which(as.character(geneExp$id) == gid ), ]

  
  xsnps <- xsnps[, grep("id", names(xsnps), invert = T)]
  xsnps <- data.frame(colnames(xsnps), t(xsnps))
  colnames(xsnps) <- c("sample", "genotype")
  
  
  
  xexp <- xexp[, grep("id", names(xexp), invert = T)]
  xexp <- data.frame(colnames(xexp), t(xexp))
  colnames(xexp) <- c("sample", "expression")
  
  df1 <- merge(xsnps, xexp, by.x = "sample", by.y = "sample")
  df1$genotype = as.factor(df1$genotype)
  mdf <- merge(df1, covar2, by.x = "sample", by.y = "sample")
  mdf$population = as.factor(mdf$population)
  mdf$gender = as.factor(mdf$gender)
  mdf$group = as.factor(mdf$group)
  
  
  mdf$presence = as.factor(mdf$genotype)
  levels(mdf$presence) = c("TE.absent", "TE.present", "TE.present")
  
  return(mdf)
}



ploteqtl <- function(x, by = "group"){
  mdf = checkExp(x)
  
  sid = as.character(x$snps)
  gid = as.character(x$gene)  
  gname = as.character(x$name)
  pval = format(as.numeric(x$pvalue), digits = 3, scientific = T)
  chisq.p = format(as.numeric(x$chisq.p), digits = 3, scientific = T)
  
  
  if(by == "genotype"){
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = genotype)) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw()
  }else if(by == "population"){
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = population)) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + scale_color_manual(values = c(rep("red", 3), "steelblue", "red"))
  }else if(by == "gender"){
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = gender)) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  }else if(by == "group"){
    g <- ggplot(mdf, aes(genotype, expression)) + geom_point(size = 3.5, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = factor(group))) + geom_boxplot(alpha = 1/2, outlier.size = NA) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + scale_color_manual(values = c("red", "blue"))
    print(g)
    slices <- table(mdf$genotype, mdf$group)
    par(mfrow=c(1,3))
    for(i in 1:nrow(slices)){
      pie(slices[i, ], col = c("red", "blue"), labels = "", border = NA)
    }   
  }else if(by == "none"){
    ggplot(mdf, aes(genotype, expression)) + geom_point(size = 3.5, alpha = 1/2, position = position_jitter(width = 0.2), col = "steelblue") + geom_boxplot(alpha = 1/2, outlier.size = NA) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
  }else{ print("No such feature for individuals") }
}


# ens2hgnc <- function(x) {
#   require(biomaRt)
#   res1 <- getMapping(x, "ensembl_gene_id", attributes = "hgnc_symbol")
# }
# 
# getMapping <- function(ids, id.type, attributes){
#   ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#   res <- getBM( attributes = c(paste(id.type), paste(attributes)), filters = id.type, values = ids, mart = ensembl )
# }


extractTFTG <- function(data){
  require(stringr)
  tf = str_extract(data, "bindingFactor=.*?\\s")
  tf = gsub("\\s$", "", tf)
  tf = gsub("bindingFactor=", "", tf)
  tf = gsub("-isoform[0-9]", "", tf)
  
  tg = str_extract(data, "hgnc=.*?;")
  tg = gsub(";$", "", tg)
  tg = gsub("hgnc=", "", tg)
  
  res = data.frame(tf = tf, target = tg)
}

extractTFTGchip <- function(data){
  require(stringr)
  tf = str_extract(data, "bindingFactor=.*?;")
  tf = gsub(";$", "", tf)
  tf = gsub("bindingFactor=", "", tf)
  tf = gsub("-isoform[0-9]", "", tf)
  
  tg = str_extract(data, "hgnc=.*?;")
  tg = gsub(";$", "", tg)
  tg = gsub("hgnc=", "", tg)
  
  res = data.frame(tf = tf, target = tg)
}



plotMDS <- function(fit, col, labels) {
  x <- fit$points[,1]
  y <- fit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="", pch = 19, col = col, bty = "n")
  legend("topright", labels, pch = 19, col=levels(as.factor(col)), cex=1)
}

fitMDS <- function(x) {
  x <- t(as.matrix(x))
  storage.mode(x) <- "integer"
  d <- dist(x)
  fit <- cmdscale(d, eig = TRUE, k = 2)
  return(fit)
}


extractExpression <- function(x, use.mean = T, show.group = T) {
  te = x$snps[1]
  gene.info = x[, c("name", "gene")]
  colnames(gene.info) = c("name", "id")
  
  expressions = as.matrix(valCols(geneExp[which(geneExp$id %in% gene.info$id),]))
  
  genotypes = valCols(snps[snps$id == te, ])
  
  indiv = unlist(genotypes)
  type0 = names(indiv)[grep(0, indiv)]
  type1 = names(indiv)[grep(1, indiv)]
  type2 = names(indiv)[grep(2, indiv)]
  
  if(use.mean){
    exp0 = rowMeans(expressions[, type0])
    exp1 = rowMeans(expressions[, type1])
    if(length(type2) <= 1){
      exp2 = expressions[, type2]
    }else{
      exp2 = rowMeans(expressions[, type2])
    }
  }else{
    exp0 = expressions[, type0]
    colnames(exp0) = paste0(colnames(exp0), "_____0")
    exp1 = expressions[, type1]
    colnames(exp1) = paste0(colnames(exp1), "_____1")
    exp2 = expressions[, type2]
    if(length(type2) > 1){
      colnames(exp2) = paste0(colnames(exp2), "_____2")
    }
  }
  
  
  if(length(type2)>0) {
    res = data.frame(exp0, exp1, exp2)
  } else {
    res = data.frame(exp0, exp1)
  }
  
  #Convert ensembl ID to hgnc Symbol
  id.map = id.map.u
  rownames(id.map) = id.map$gene
  id.map = id.map[rownames(res), ]
  rownames(res) = id.map$name
  
  #population group pie chart
  if(show.group){
    mdf = checkExp(x[1,])
    slices <- table(mdf$genotype, mdf$group)
    par(mfrow=c(1,3))
    for(i in 1:nrow(slices)){
      pie(slices[i, ], col = c("red", "blue"), labels = "", border = NA)
    }
    mtext(te)
  }
    
  
  return(res)
}


extractTarget <- function(x, isOriginal = F) {
  te = x$snps[1]
  if(isOriginal){
    gene.info = x[, c("name", "gene")]
  } else {
    gene.info = x[, c("target", "gene")]
  }
  colnames(gene.info) = c("name", "id")
  expressions = as.matrix(valCols(geneExp[which(geneExp$id %in% gene.info$id),]))
 
  genotypes = valCols(snps[snps$id == te, ])
  
  indiv = unlist(genotypes)
  type0 = names(indiv)[grep(0, indiv)]
  type1 = names(indiv)[grep(1, indiv)]
  type2 = names(indiv)[grep(2, indiv)]
  
  
  exp0 = rowMeans(expressions[, type0])
  exp1 = rowMeans(expressions[, type1])
  exp2 = rowMeans(expressions[, type2])
  
  if(length(type2)>0) {
    res = data.frame(exp0, exp1, exp2)
  } else {
    res = data.frame(exp0, exp1)
  }
  
  #Convert ensembl ID to hgnc Symbol
  id.map = id.map.u
  rownames(id.map) = id.map$gene
  id.map = id.map[rownames(res), ]
  rownames(res) = id.map$name
  
  return(res)
}


zscale <- function(x){
  x = as.matrix(x)
  z = (x - rowMeans(x))/apply(x, 1, sd)
}

medScale <- function(x){
  x = as.matrix(x)
  z = x - rowMeans(x)
}

myHeatmap <- function(data, 
                      scaled = T, 
                      title = "", 
                      hmcol = colorRampPalette(c("green", "darkgreen", "black","darkred", "red"))(100), 
                      labCol = NULL
                      ){
  require(gplots)
  
  if(scaled){
    if(ncol(data) > 2){
      data = zscale(data)
    }else{
      #data = medScale(data)
      data = medScale(data)
    }
  }else{
    data = as.matrix(data)
  }
  
  heatmap.2(data, 
            col = hmcol,
            Colv = FALSE,
            trace = "none", 
            scale = "none",
            labRow = rownames(data),
            labCol = labCol,
            margin = c(10, 6),
            keysize = 1.5,
            density.info = "none",
            key.title = "z-score",
            key.xlab = "z-score",
            main = title)
}
    