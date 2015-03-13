
getChisq <- function(x, verbose = F){
  mdf = checkExp(x, verbose)
  tb = table(mdf$group, mdf$presence)
  chisq.test(tb)$p.value
}

ploteqtl <- function(x, by = "genotype" ){
  mdf = checkExp(x)
  
  sid = as.character(x$snps)
  gid = as.character(x$gene)  
  pval = format(as.numeric(x$pvalue), digits = 3, scientific = T)
  chisq.p = format(as.numeric(x$chisq.p), digits = 3, scientific = T)
  
    
  if(by == "genotype"){
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = genotype)) + labs(title = paste(sid, "vs", gid)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw()
  }else if(by == "population"){
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = population)) + labs(title = paste(sid, "vs", gid)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw() + scale_color_manual(values = c(rep("red", 3), "steelblue", "red"))
  }else if(by == "gender"){
    ggplot(mdf, aes(genotype, expression)) + geom_boxplot() + geom_point(size = 4, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = gender)) + labs(title = paste(sid, "vs", gid)) + xlab(paste("P=", pval, " Chisq.P=", chisq.p)) + theme_bw()
  }else{ print("No such feature for individuals") }
}

checkExp <- function(x, verbose = F){
  #lapply(1:5, function(x) ploteqtl(alleqtl[x, ]))
  #color points by "genotype" or "population"
  if(verbose) print(x)
  
  require(ggplot2)
  sid = as.character(x$snps)
  gid = as.character(x$gene)  
  pval = as.numeric(x$pvalue)
  
  xsnps <- snps[which(as.character(snps$id) == sid ), ]
  xexp <- geneExp[which(as.character(geneExp$id) == gid ), ]
  
  xsnps <- xsnps[, grep("id", names(xsnps), invert = T)]
  xsnps <- data.frame(colnames(xsnps), t(xsnps))
  colnames(xsnps) <- c("sample", "genotype")

  
  
  xexp <- xexp[, grep("id", names(xexp), invert = T)]
  xexp <- data.frame(colnames(xexp), t(xexp))
  colnames(xexp) <- c("sample", "expression")
  
  df1 <- merge(xsnps, xexp, by.x = "sample", by.y = "sample")
  df1$genotype = as.factor(df1$genotype)
  #df2 = data.frame(sample = colnames(covar2)[grep("id", colnames(covar), invert = T)], t(covar[, grep("id", colnames(covar), invert = T)]))
  #if(ncol(df2)==2) 
  #  colnames(df2) = c("sample", "gender")
  #else
  #  colnames(df2) = c("sample", "population", "gender")
  mdf <- merge(df1, covar2, by.x = "sample", by.y = "sample")
  #if(ncol(df2) > 2) mdf$population = as.factor(mdf$population)
  mdf$population = as.factor(mdf$population)
  mdf$gender = as.factor(mdf$gender)

  
  mdf$presence = as.factor(mdf$genotype)
  levels(mdf$presence) = c("TE.absent", "TE.present", "TE.present")
  
  return(mdf)
}




runPop <- function(x, af.cut, useModel = modelLINEAR, par.plot = 100, verbose = FALSE, output_file_name = tempfile(), base.dir = ".", pvOutputThreshold = 1e-2)
{
	#path to input
	SNP_file_name = paste0("./temp3/snps", x, ".txt")
	snps_location_file_name = "./input/snploc.txt"

	expression_file_name = paste0("./temp3/geneExp", x, ".txt")
	gene_location_file_name = paste0("./temp3/geneloc", x, ".txt")
	covariates_file_name = paste0("./temp3/covariates", x, ".txt")

	# read some input
	snps = read.table(SNP_file_name, header = T, stringsAsFactors = F)
	snploc = read.table(snps_location_file_name, header = T, stringsAsFactors = F)
	
	require(dplyr)
	# getAF
	msnps <- as.matrix(snps[, grep("id", names(snps), invert = T)])
	rownames(msnps) <- as.character(snps[, grep("id", names(snps))])
	res <- apply(msnps, 1, table)
	res1 <- sapply(res, getAF)
	af.df <- data.frame(names(res1), res1)
	names(af.df) = c("id", "af")
	af.df$id = as.character(af.df$id)
	
	# AF filter
	keep = which(af.df$af > af.cut)
	res = snps[keep, ]
	res1 = snploc[keep, ]
	write.table(res, "./temp/temp.snps", quote = F, row.names = F, col.names = F, sep = "\t")
	write.table(res1, "./temp/temp.snploc", quote = F, row.names = F, col.names = F, sep = "\t")
	
	# Filtered input
	SNP_file_name = "./temp/temp.snps"
	snps_location_file_name = "./temp/temp.snploc"
	
	# Get results + histogram
	meh <- runME( SNP_file_name, snps_location_file_name, expression_file_name, gene_location_file_name, covariates_file_name, useModel = useModel, verbose = verbose)
	
	# Get results + qqplot
	meq <- runME( SNP_file_name, snps_location_file_name, expression_file_name, gene_location_file_name, covariates_file_name, useModel = modelLINEAR, "qqplot", verbose = F)
	
	res <- list(meh = meh, meq = meq)
}

getAF <- function(x){
    x <- as.vector(x)
    #print(x)
    #print(sum(x))
    if(length(x) == 1)
    {
        res = 0
    }
    else if( length(x) == 2 )
    {
        res <- x[2]/(sum(x)*2)
    }
    else
    {
        res <- (x[2] + 2*x[3])/(sum(x)*2)
    }
    return(res)
}

runME <- function( SNP_file_name, snps_location_file_name, expression_file_name, gene_location_file_name, covariates_file_name, useModel = modelLINEAR, par.plot = 100, verbose = FALSE, output_file_name = tempfile(), base.dir = ".", pvOutputThreshold = 1e-2)
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
	noFDRsaveMemory = FALSE);

	unlink(output_file_name);
	
	return(me)
}
