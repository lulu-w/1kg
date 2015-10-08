loadPvals <- function(pvalfile){
  require(data.table)
  pvals = fread(pvalfile, header = F)
  setnames(pvals, c("snp", "te", "statistic", "pval", "FDR", "beta", "chr", "pos"))
  pvals = as.data.frame(pvals)
  pvals = pvals[, c("snp", "pval", "chr", "pos", "te")]
  pvals = pvals[order(pvals$chr, pvals$pos), ]
}

myQQplot = function(pvalfile, permutedfile, n_probs=100000, p_cut=0.01, col="black", title=""){
  require(data.table)
  pvals = fread(pvalfile, header = F)
  ppvals = fread(permutedfile, header = F)
  setnames(pvals, c("snp", "te", "statistic", "pval", "FDR", "beta"))
  setnames(ppvals, c("snp", "te", "statistic", "pval", "FDR", "beta"))
  observed = pvals$pval
  pobserved = ppvals$pval
  med_observed=median(observed)
  med_pobserved=median(pobserved)
  observed = -log10(observed)
  pobserved = -log10(pobserved)
  pop_name=title
  ptitle=paste(title, "permuted")
  title=paste(title, "observed")

  probs=1:n_probs/n_probs
  obs_q=quantile(observed, probs=probs)
  pobs_q=quantile(pobserved, probs=probs)
  exp_q=qunif(probs, max=p_cut)
  med_expected=median(exp_q)
  exp_q=sort(-log10(exp_q))
  lambda=med_observed/med_expected
  plambda=med_pobserved/med_expected
  
  min_exp_q=round(exp_q[1])
  max_exp_q=round(exp_q[length(exp_q)])
  min_obs_q=round(obs_q[1])
  max_obs_q=round(obs_q[length(obs_q)])
  
  par(cex=0.7)
  tick_space=10
  plot(exp_q, obs_q, col=col, xlab=expression(-log[10]*(Expected)), ylab=expression(-log[10]*(Observed)), main=pop_name, ylim=c(min_obs_q, (max_obs_q %/% tick_space + 1)*tick_space), xlim=c(min_exp_q, max_exp_q+1), axes=F)
  points(exp_q, pobs_q, col="darkgrey")
  
  legend("top", 
         col=c(col, "darkgrey"), 
         pch=1, 
         legend=c(as.expression(bquote(.(title) * ", " * lambda == .(round(lambda, 3)))), as.expression(bquote(.(ptitle) * ", " * lambda == .(round(plambda, 3))))), 
         bty='n')
  
  axis(side=1, pos=min_exp_q)
  axis(side=2, pos=min_obs_q, at=c(min_obs_q, seq((min_obs_q %/% tick_space + 1)*tick_space, (max_obs_q %/% tick_space + 1)*tick_space, by=tick_space)))
  abline(0,1)
  abline(h=8, col="black", lty=2)
}


plotManhattan <- function(pvals, col=c("black", "red"), title="", highlight=NA, by="none"){  
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
  #pvals$chr = factor(pvals$chr, levels = paste0("chr", 1:chrNum))
  pvals$chr = factor(pvals$chr, labels = paste0("chr", 1:chrNum))
  
  require(ggplot2)
  
  if(is.na(highlight[1])){
      if(by=="none"){
        ggplot(pvals, aes(x = pos, y = -log10(pval))) + 
          geom_point(aes(color = chr)) +
          scale_color_manual(values = rep(col, chrNum/2), guide = F) + 
          theme_bw() + 
          geom_text(data = data.frame(text = 1:chrNum, x = xlab.pos), aes(x = x, y = 0, label = text), size = 3)+
          theme(panel.grid.major=element_blank(), 
                panel.grid.minor=element_blank(), 
                axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) +
          labs(title = title) + 
          xlab("chromosomes") + 
          ylab(bquote(-log[10]*italic(P)*"-value")) +
          coord_cartesian(ylim=c(-2, 70)) +
          geom_hline(aes(yintercept=8), color=col2inv(col[2]))
      }else{
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
  }else{
      if(length(col)==2){ col = c(col, "green")}
      pvals$group = ifelse(pvals$snp %in% highlight, "z", pvals$group)
      ggplot(pvals, aes(x = pos, y = -log10(pval))) + 
        geom_point(aes(color = group, shape = te, size = -log10(pval))) +
        scale_shape_manual(values = c(16, 15, 17, 18)) +
        # scale_color_manual(values = rep(col, chrNum/2), guide = F) + 
        scale_color_manual(values = col, guide = F) + 
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
  
}



histLogFreq = function(d, breaks=1000, col="black", title=""){
  hist.data = hist(d, breaks=breaks, plot=F)
  hist.data$counts = log10(hist.data$counts)
  plot(hist.data, ylab=bquote(-log[10]*"Frequency"), ylim=c(0,5), xlab=bquote(-log[10]*italic(P)*"-value"), border=NA, col=col, main=title, xlim=c(0, 60))
  abline(v=8, col="black", lty=2)
}


addGTEx <- function(filename){
    data = read.delim(filename, sep = "\t", header = T )
    data = data[,c(1:6,9)]
    colnames(data)[1] = "snp"
    return(data)
}


# regressRecomb <- function(pvals, recomb){
#     recomb = recomb[, 1:5]
#     require(dplyr)
#     require(ggplot2)
#     bin_size = 1e6
#     pvals$bin_start = as.integer((pvals$pos %/% bin_size) * bin_size)
#     pvals$bin_end = as.integer(pvals$bin_start + bin_size)
#     pvals$window = paste0(pvals$chr, "_", pvals$bin_start, "_", pvals$bin_end)
#     recomb$window = paste0(recomb[, 1], "_", recomb$chromStart, "_", recomb$chromEnd)
#     pvals = left_join(pvals, recomb, by = "window")
#     g <- ggplot(pvals, aes(x=decodeAvg)) + geom_histogram(binwidth = 0.05)
#     print(g)
#     # Replace NA recomb rate with -1
#     # pvals$decodeAvg = ifelse(is.na(pvals$decodeAvg), -1, pvals$decodeAvg)
#     # ggplot(pvals, aes(x=decodeAvg, y = pval)) + geom_point()
#     #g <- ggplot(pvals, aes(x=decodeAvg, y = -log10(pval))) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~ te, nrow=2)
#     #print(g)
#     df = pvals[complete.cases(pvals), ]
#     #par(mfrow=c(2, 2))
#     plotfit(df[which(df$te=="all"), "decodeAvg"], df[which(df$te=="all"), "pval"], title="all")
#     plotfit(df[which(df$te=="Alu"), "decodeAvg"], df[which(df$te=="Alu"), "pval"], title="Alu")
#     plotfit(df[which(df$te=="L1"), "decodeAvg"], df[which(df$te=="L1"), "pval"], title="L1")
#     plotfit(df[which(df$te=="SVA"), "decodeAvg"], df[which(df$te=="SVA"), "pval"], title="SVA")
# }


loadDecode <- function(path, chr="all"){
    if(chr=="all"){
        recomb = lapply(1:22, function(i) {
            filename = paste0(path, "decodeSexAveraged_Recomb_hg19_chr", i, "_20150714.txt")
            res = read.delim(filename, sep = "\t", header = F, comment.char = "#", skip=1)
        })
        recomb = do.call(rbind, recomb)
    }else{
        filename = paste0(path, "decodeSexAveraged_Recomb_hg19_chr", chr, "_20150714.txt")
        recomb = read.delim(filename, sep = "\t", header = F, comment.char = "#", skip=1)
    }
    return(recomb)
}


pvalJoinRecomb <- function(pvals, recomb){
    require(dplyr)
    require(ggplot2)
    bin_size = 1e4
#     pvals$bin_start = as.integer((pvals$pos %/% bin_size) * bin_size)
#     pvals$bin_end = as.integer(pvals$bin_start + bin_size)
#     pvals$window = paste0(pvals$chr, "_", pvals$bin_start, "_", pvals$bin_end)

    pvals_list = split(pvals, f=pvals$chr)
    recomb_list = split(recomb, f=recomb$chr)
    name_list=names(pvals_list)

    res = mapply(function(pvals, recomb, idx){
        recomb_max = recomb$chromEnd[nrow(recomb)]
        pvals_max = pvals$pos[nrow(pvals)]
        if(pvals_max > recomb_max){
            cutpoints = as.integer(c(0, recomb$chromStart, recomb_max, pvals_max))
        }else{
            cutpoints = c(0, recomb$chromStart, recomb_max)
        }
        pvals$posBin=cut(pvals$pos, cutpoints)
        recomb$posBin=cut(recomb$chromEnd, cutpoints)
        pvals=left_join(pvals, recomb, by="posBin")
        return(pvals)
    }, pvals_list, recomb_list, 1:length(name_list), SIMPLIFY=F)

    res=do.call(rbind, res)
    return(res)
}


plotRecomb <- function(pvals, title="", binwidth=1, yscale=NA){
    if(!is.na(yscale) & yscale=="log2"){ pvals$decodeAvg = log2(pvals$decodeAvg)}
    g = ggplot(pvals, aes(x=decodeAvg)) + 
      geom_histogram(binwidth = binwidth) + 
      labs(title=title) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank())
    if(!is.na(yscale) & yscale=="log2"){
      g = g+labs(x="log2 decodeAvg")
    }
    print(g)
}

plotRegress <- function(pvals, title="", col=rgb(0,0,1,0.2)){
    df = pvals[complete.cases(pvals), ]
    par(mfrow=c(2,2))
    plotfit(df[which(df$te=="all"), "decodeAvg"], df[which(df$te=="all"), "pval"], title=paste("all from", title), col = col)
    plotfit(df[which(df$te=="Alu"), "decodeAvg"], df[which(df$te=="Alu"), "pval"], title=paste("Alu from", title), col = col)
    plotfit(df[which(df$te=="L1"), "decodeAvg"], df[which(df$te=="L1"), "pval"], title=paste("L1 from", title), col = col)
    plotfit(df[which(df$te=="SVA"), "decodeAvg"], df[which(df$te=="SVA"), "pval"], title=paste("SVA from", title), col = col)
}


plotfit <- function(x, y, title = "", yscale = "-log10", col = col){
    par(cex = 0.7)
    if(yscale=="-log10") { y = -log10(y)}
    m = lm(y ~ x)
    # ~ space, * juxtapose, .() use variable content
    eq = bquote(atop(italic(y) == .(round(coef(m)[1], 3)) + ~"("*.(round(coef(m)[2], 3))*")" %.% italic(x), italic(r)^2 == .(round(summary(m)$r.squared, 4))))
    
    plot(x, y, pch = 16, col = col, xlab = "decodeAvg", ylab = "-log10(pval)", main = title)
    abline(m, col = col2inv(col))
    mtext(eq, side = 3, line = -3, cex = 0.7)
    #legend("topright", (eq))
}

col2alpha <- function(col, alpha){
    res = lapply(col, function(x){
        rgb(col2rgb(x)[1,1]/255, col2rgb(x)[2,1]/255, col2rgb(x)[3,1]/255, alpha)
    })
    res = unlist(res)
}

col2inv <- function(col){
    res = lapply(col, function(x){
      rgb((255-col2rgb(x)[1,1])/255, (255-col2rgb(x)[2,1])/255, (255-col2rgb(x)[3,1])/255)
    })
    res = unlist(res)
}


ggplotRegress <- function(pvals, title="", col=rgb(0,0,1,0.2), bin.data=NA){
    pvals=pvals[complete.cases(pvals),]
    pvals$pval=-log10(pvals$pval)
    if(is.na(bin.data)){
        ggplot(pvals, aes(x=decodeAvg, y=pval)) + 
          geom_point(col=col, alpha=0.2) +
          geom_smooth(method="lm", col = col2inv(col)) +
          facet_wrap(~te, nrow=2) + 
          theme_bw() +
          theme(panel.grid.major = element_blank())
    }else if(bin.data=="decodeAvg"){
        cutpoints = c(0, 5, 10, max(pvals$decodeAvg))
        pvals$decodeAvgBin = cut(pvals$decodeAvg, cutpoints, include.lowest = T)
        ggplot(pvals, aes(x=decodeAvg, y=pval)) + 
          geom_point(col=col, alpha=0.2) +
          geom_smooth(method="lm", col = col2inv(col)) +
          facet_wrap(decodeAvgBin~te, ncol=4, scale = "free") + 
          theme_bw() +
          theme(panel.grid.major = element_blank())
    }
}


addBins <- function(data, by, bin_size){
    data$bin_start=as.integer((data[, by] %/% bin_size) * bin_size)
    data$bin_end=as.integer(data$bin_start + bin_size)
    return(data)
}


pvalJoinGaps <- function(pvals, gaps){
  require(dplyr)
  require(ggplot2)
  
  pvals_list = split(pvals, f=pvals$chr)
  name_list=names(pvals_list)
  gaps = gaps[which(gaps$chrom %in% name_list), ]
  gaps_list = split(gaps, f=gaps$chrom)
  
  
  res = mapply(function(pvals, gaps, idx){
    cutpoints=sort(c(gaps$chromStart, gaps$chromEnd))
    pvals$posBin=cut(pvals$pos, cutpoints)
    keep=levels(pvals$posBin)[c(1,3,5)]
    pvals=pvals[which(pvals$posBin %in% keep), ]
    
    gaps$posBin=cut(gaps$chromEnd, cutpoints)
    pvals=left_join(pvals, gaps, by="posBin")
    return(pvals)
  }, pvals_list, gaps_list, 1:length(name_list), SIMPLIFY=F)
  
  res=do.call(rbind, res)
  return(res)
}


summaryMappability <- function(pvals, feature_file_path){
    pvals_list=split(pvals, f=pvals$chr)
    res=mapply(function(pval, chr_name){
        require(dplyr)
        require(data.table)
        cur_feat_file=paste0(feature_file_path, "wgEncodeDukeMapabilityUniqueness20bp_", chr_name, ".bedGraph")
        #feat=read.delim(cur_feat_file, header=F, sep="\t")
        feat=fread(cur_feat_file, sep="\t", header = F)
        setnames(feat, c("chr", "chromStart", "chromEnd", "uniqueness"))
        # Assuming ALL data are sorted by position
        pval_max=pval$pos[length(pval$pos)]
        feat_max=feat$chromEnd[length(feat$chromEnd)]
        pval_min=pval$pos[1]
        feat_min=feat$chromStart[1]
        minSmall=min(pval_min, feat_min)
        if(pval_max > feat_max){
          if(minSmall==0){
              cutpoints = as.integer(c(feat$chromStart, feat_max, pval_max))
          }else{
              cutpoints = as.integer(c(0, feat$chromStart, feat_max, pval_max))
          }
        }else{
          if(minSmall==0){
              cutpoints = as.integer(c(feat$chromStart, feat_max))
          }else{
              cutpoints = as.integer(c(0, feat$chromStart, feat_max))
          }
        }
        pval=as.data.table(pval)
        pval[, posBin:=cut(pos, cutpoints)]
        setkey(pval, posBin)
        feat[, posBin:=cut(feat$chromEnd, cutpoints)]
        setkey(feat, posBin)
        pval=left_join(pval, feat, by="posBin")
        tb=table(pval$uniqueness)
        fg_freq=as.vector(tb)
        names(fg_freq)=names(tb)
        feat[, bases:=chromEnd-chromStart]
        # Remove the first and last chunk of bases in feature track
        feat$bases[1]=0
        feat$bases[nrow(feat)]=0
        feat[, sumBases:=sum(bases), by=uniqueness]
        bg_freq=unique(feat$sumBases)
        bg_freq_names=unique(feat$uniqueness)
        names(bg_freq)=bg_freq_names
        bg_freq=bg_freq[order(bg_freq_names)]
        freq=rbind(fg_freq, bg_freq)
        res=list(freq, pval)
        return(res)
    }, pvals_list, names(pvals_list), SIMPLIFY=F)
    #res=do.call(rbind, res)
    return(res)
}


# pvalJoinMappability <- function(pvals, feature_file_path){
#     pvals_list=split(pvals, f=pvals$chr)
#     res=mapply(function(pval, chr_name){
#         require(dplyr)
#         cur_feat_file=paste0(feature_file_path, "wgEncodeDukeMapabilityUniqueness20bp_", chr_name, ".bedGraph")
#         feat=read.delim(cur_feat_file, header=F, sep="\t")
#         colnames(feat)=c("chr", "chromStart", "chromEnd", "uniqueness")
#         # Assuming ALL data are sorted by position
#         pval_max=pval$pos[length(pval$pos)]
#         feat_max=feat$chromEnd[length(feat$chromEnd)]
#         pval_min=pval$pos[1]
#         feat_min=feat$chromStart[1]
#         minSmall=min(pval_min, feat_min)
#         if(pval_max > feat_max){
#           if(minSmall==0){
#               cutpoints = as.integer(c(feat$chromStart, feat_max, pval_max))
#           }else{
#               cutpoints = as.integer(c(0, feat$chromStart, feat_max, pval_max))
#           }
#         }else{
#           if(minSmall==0){
#               cutpoints = as.integer(c(feat$chromStart, feat_max))
#           }else{
#               cutpoints = as.integer(c(0, feat$chromStart, feat_max))
#           }
#         }
#         pval$posBin=cut(pval$pos, cutpoints)
#         feat$posBin=cut(feat$chromEnd, cutpoints)
#         pval=left_join(pval, feat, by="posBin")
#         return(pvals)
#     }, pvals_list, names(pvals_list), SIMPLIFY=F)
#     res=do.call(rbind, res)
#     return(res)
# }


plotUniqueness <- function(data, col=NA, title=""){
    pvals=data[[2]]
    count=data[[1]]
    if(is.na(col)){
        hmcol=heat.colors(12)
    }else{
        pal=colorRampPalette(c("grey", col))
        hmcol=pal(12)
    }
  
    pvals$logPval=-log10(pvals$pval)
    pvals=pvals[order(pvals$logPval), ]
    
    uniquenessVals=names(table(pvals$uniqueness, useNA="ifany"))
    uniqueness_m=matrix(rep(uniquenessVals, each=nrow(pvals)), nrow=nrow(pvals))
    colnames(uniqueness_m)=uniquenessVals
    uniqueness_m=apply(uniqueness_m, 2, function(x){ x==pvals$uniqueness})*1
    
    par(mfrow=c(7,1), mar=c(1,4,1,1), oma=c(2,0,2,0), cex=0.7)
    barplot(count[1, ]/count[2, ], col=col, yaxt='n', ylab="ratio pval/bg", names.arg=paste(c(1:4, ">4"), "time(s)"))
    
    for(i in 1:ncol(uniqueness_m)){
      mapTimes=ifelse(colnames(uniqueness_m)[i]==0, ">4", round(1/as.numeric(colnames(uniqueness_m)[i])))
      image(as.matrix(uniqueness_m[, i]), ylab=paste(mapTimes, "time(s)"), axes=F, col=hmcol)
    }
    plot(1:nrow(pvals), pvals$logPval, type="n", ylab="ordered p-value", xaxt="n")
    lines(1:nrow(pvals), pvals$logPval, col=col)   
    mtext("ordered -log10 p-value", side=1, cex=0.8)
    mtext(title, side=3, cex=0.8, outer=T)
}
