source("functions-analysis.R")

# test_rs = data.frame(snp=c("rs77030489", "rs79709321"), chr=c("chr16", "chr21"), pos=c(46423036, 11128387), te=c("SVA", "L1"), pval=c("", ""))
# test_rs = addBins(test_rs, by="pos", bin_size=1e6)
# pvals=test_rs

pvals = loadPvals("./analysis_input/chrAll_AFR_pval.txt")
pvals = pvals[order(pvals$pval), ]
top_n=2
pvals=pvals[1:top_n, ]
pvals = addBins(pvals, by="pos", bin_size=1e6)

population="AFR"

####
# Expression data
exprs=read.table("./analysis_input/presenceExp.txt", header = T)
rownames(exprs)=exprs$id
exprs = exprs[,-1]

####
# Covariates data
covar=read.table("./abilbeast_input/sampleInfo.txt", header = T)

####
# SNP genotype data
genotype_files=paste0("./abilbeast_input/", gsub("chr", "", pvals$chr), "_", (pvals$bin_start+1), "_", pvals$bin_end, ".genotype")

pvals=split(pvals, f=1:nrow(pvals))
res=mapply(function(pval, genotype){
    require(dplyr)
    require(ggplot2)
    
    te=unlist(pval["te"])
    snp=unlist(pval["snp"])
    pvalue=unlist(pval["pval"])
    chr=unlist(pval["chr"])
    pos=unlist(pval["pos"])
    
    ex=unlist(exprs[paste(te), ])
    df=data.frame(id=names(ex), expression=ex)
    
    genotype_data=read.delim(genotype, header = T, sep = "\t")
    genotype_data=genotype_data[which(genotype_data[, 1]==snp), ]
    genotype_data=unlist(genotype_data[, -1])
    df$genotype=genotype_data
    
    mdf=left_join(df, covar, by="id")
    mdf=mdf[which(mdf$group==population), ]
    
    g=ggplot(mdf, aes(x=factor(genotype), y=expression)) + 
      geom_point(size = 2, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = factor(genotype))) +
      geom_boxplot(alpha = 1/2, outlier.size = 0) + 
      labs(title = paste(chr, ":", pos, "\n\n", snp, "vs", te)) + 
      xlab(paste("SNP genotype\n\nP.value = ", format(pvalue, digits=4))) + 
      ylab(paste(te, "insertion counts")) +
      theme_bw() + 
      theme(panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank(),
            plot.title=element_text(size=rel(0.9)))
    return(g)
  
}, pvals, genotype_files, SIMPLIFY=F)

res=lapply(res, print)
