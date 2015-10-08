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


plotManhattan <- function(pvals, col=c("black", "red"), title="", highlight=NA, by="none"){  
  if(length(grep("chr", pvals$chr)) == 0){
      pvals$chr = paste0("chr", pvals$chr)
  }
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
  pvals$chr = factor(pvals$chr, labels = paste0("chr", 1:chrNum))
  #pvals = pvals[order(pvals$chr, pvals$pos), ]
  
  require(ggplot2)
  
  if(is.na(highlight[1])){
    if(by=="none"){
      ggplot(pvals, aes(x = pos, y = -log10(pval))) + 
        geom_point(aes(color = group)) +
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