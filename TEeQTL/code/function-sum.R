

checkExp <- function(x, verbose = F){
  #lapply(1:5, function(x) ploteqtl(alleqtl[x, ]))
  #color points by "genotype" or "population"
  if(verbose) print(x)
  
  require(ggplot2)
  sid = as.character(x$snps)
  gid = as.character(x$gene)  
  pval = as.numeric(x$pvalue)
  
  xsnps <- snps[which(snps$ID == sid ), ]
  xexp <- geneExp[which(as.character(geneExp$id) == gid ), ]
  
  
  xsnps <- xsnps[, grep("ID", names(xsnps), invert = T)]
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
  #gname = as.character(x$name)
  pval = format(as.numeric(x$pvalue), digits = 3, scientific = T)
  #chisq.p = format(as.numeric(x$chisq.p), digits = 3, scientific = T)
  
  ggplot(mdf, aes(genotype, expression)) + geom_point(size = 3.5, alpha = 1/2, position = position_jitter(width = 0.2), aes(col = factor(group))) + geom_boxplot(alpha = 1/2, outlier.size = NA) + labs(title = paste(sid, "vs", gid)) + xlab(paste("P=", pval)) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + scale_color_manual(values = c("red", "blue"))

}
