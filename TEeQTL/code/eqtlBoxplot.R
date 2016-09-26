eqtlBoxplot <- function(d, popName, snploc, snps, exprs, color = "black", yminAdd=0, ymax=NULL){
  require(data.table)
  require(ggplot2)
  
  snp = as.character(d[1])
  geneName = as.character(d[2])
  gene = as.character(d[3])
  stat = as.numeric(d[4])
  pval = as.numeric(d[5])
  FDR = as.numeric(d[6])
  beta = as.numeric(d[7])
  
  idx = which(snploc$snp == snp) # this version of data.table does not work with which()
  dloc = unlist(snploc[idx, ])
  
  genotype = unlist(snps[which(snps[, 1] == snp), ])
  genotype = data.frame(id = names(genotype)[2:length(genotype)], allele = genotype[2:length(genotype)])
  
  expression = unlist(exprs[which(exprs$id == gene), ])
  expression = data.frame(id = names(expression)[2:length(expression)], level = expression[2:length(expression)])
  
  df = merge(genotype, expression, by.x = "id", by.y = "id")
  #df$allele = factor(df$allele, level = c(0, 1, 2), labels = c(paste0(dloc[4], dloc[4]), paste0(dloc[4], dloc[5]), paste0(dloc[5], dloc[5])))
  #df$allele = factor(df$allele)
  df$level = as.numeric(df$level)
  df$numallele = as.numeric(df$allele)
  N = length(levels(df$allele))
  #ylim = max(df$level)+ymaxAdd
  #ymin = min(df$level)+yminAdd
  
  
  p = ggplot(data = df, aes(x = allele, y = level)) +
    geom_point(aes(x = allele, y = level), position=position_jitter(0.2), col = color) +
    geom_boxplot(alpha = 0.3, outlier.size = 0) +
    labs(x = paste0("P=", format(pval, 3, scientific = T), "\n FDR=", format(FDR, 3, scientific = T), "\n slope=", round(beta, 3)), 
         y = paste(geneName, " expression level"), 
         title = paste(snp, " (", dloc[2], ":", dloc[3],")\n in", popName)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title=element_text(size = rel(0.7)),
          axis.title=element_text(size = rel(0.8)))
  if(!is.null(ymax)){
      p = p+coord_cartesian(ylim = ymax)
  }
  print(p)
}