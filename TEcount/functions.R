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


plotManhattan <- function(pvals, col=c("black", "red"), title=""){
  colnames(pvals) = c("snp", "te", "statistic", "pval", "FDR", "beta", "chr", "pos")
  pvals = pvals[, c("snp", "pval", "chr", "pos", "te")]
  pvals = pvals[order(pvals$chr, pvals$pos), ]
  
  chrNum = length(levels(factor(pvals$chr)))
  pvals$pos = round(pvals$pos/10000, 0)
  pvals$group = "odd"
  
  xlab.pos <- NULL
  
  for (i in 1:chrNum){ 
    ndx <- which(pvals$chr==paste0("chr", i))
    lstMrk <- max(pvals[ndx, "pos"])
    xlab.pos <- c(xlab.pos, lstMrk)
    if (i < chrNum) ndx2 <- which(pvals[, "chr"]==paste0("chr", i+1))
    if (i < chrNum) pvals[ndx2, "pos"] <- pvals[ndx2, "pos"] + lstMrk
    if (i %% 2 == 0) {
      pvals$group[ndx] = "even"
      pvals$col[ndx] = gsub("odd", "even", pvals$col[ndx])
    }
    
  }
  xlab.pos = 0.5*(c(0, xlab.pos[-length(xlab.pos)]) + xlab.pos)
  #Order chromosomes
  pvals$chr = factor(pvals$chr, levels = paste0("chr", 1:chrNum))
  
  require(ggplot2)
  ggplot(pvals, aes(x = pos, y = -log10(pval))) + 
    geom_point(aes(color = chr, shape = te, size = -log10(pval))) +
    scale_shape_manual(values = c(16, 15, 17, 18)) +
    scale_color_manual(values = rep(col, chrNum/2), guide = F) + 
    theme_bw() + 
    geom_text(data = data.frame(text = 1:chrNum, x = xlab.pos), aes(x = x, y = -2, label = text), size = 3)+
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = title) + 
    xlab("chromosomes")
}

getAF = function(d){
    id=d$id
    N=2*(ncol(d)-1)
    p=rowSums(d[, 2:ncol(d)])
    p=p/N
    af=data.frame(id=id, af=p)
}