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



ploteqtl <- function(x, by = "genotype" ){
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
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = group)) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) #+ scale_fill_manual()
  }else if(by == "none"){
    ggplot(mdf, aes(genotype, expression)) + geom_point(size = 3.5, alpha = 1/2, position = position_jitter(width = 0.2), col = "steelblue") + geom_boxplot(alpha = 1/4) + labs(title = paste(sid, "vs", gname)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
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


plotHeatmap <- function(data){
  require(gplots)
  #mat = as.matrix(valCols(data))
  mat = as.matrix(data)
#   log.mat = log2(mat + abs(min(mat)) + 1)
#   mat = log.mat
  
  zmat <- (mat - rowMeans(mat))/apply(mat, 1, sd)
  #require(ggplot2)
  
#   zmat1 <- zmat[, 1:2]
#   zmat2 <- zmat[, 3:4]
#   med1 <- median(zmat1)
#   med2 <- median(zmat2)
#   zmat1 = zmat1 - med1
#   zmat2 = zmat2 - med2
#   zmat = cbind(zmat1, zmat2)
  
  heatmap.2(zmat, 
            col = hmcol,
            #Rowv = FALSE,
            Colv = FALSE,
            trace = "none", 
            scale = "none",
            #labRow = data[, grep("Symbol", colnames(data))],
            labRow = "",
            margin = c(15, 6),
            keysize = 1.5,
            density.info = "none",
            key.title = "z-score",
            key.xlab = "z-score",
            main = "")
}

extractTarget <- function(x, te) {
  gene.info = data.frame(id = x$gene, symbols = x$target)

  expressions = valCols(geneExp[which(geneExp$id %in% gene.info$id),])
 
  genotypes = valCols(snps[snps$id == te, ])
  
  indiv = unlist(genotypes)
  indiv = colnames(genotypes)[order(indiv)]
  expressions = expressions[, indiv] 
  
  return(expressions)
}