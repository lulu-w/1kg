getGenoPhenoMatrix = function(pop_sel, te_type){
    pop_sel_id = sampleInfo$id[which(sampleInfo$population %in% pop_sel)]
    pop_sel_id = sort(pop_sel_id)
    
    # selected TE counts and phenotype Euclidean distance
    counts = unlist(teCount[which(teCount$id == te_type), pop_sel_id])
    delta_y = as.matrix(dist(as.matrix(counts), diag = T))
    
    # selected genetic distances
    sel = match(pop_sel_id, names(grm)) # want to sort names(grm) as pop_sel_id, i.e. find id in pop_sel_id one by one in names(grm), so use it as value to match against
    if(sum(is.na(sel)) > 0){
        print("Error!! Could not find individual in genetic relationship matrix")
        return()
    }
    sel_grm = as.matrix(grm[sel, sel, with=F])
    res = list(gdist = sel_grm, pdist = delta_y)
}

getMatrixColors = function(gdist, col = c("blue", "red")){
    pop_sel_id = colnames(gdist)
    pop_sel_pop = sampleInfo$pop[match(pop_sel_id, sampleInfo$id)]
    pairs_pop = expand.grid(pos1 = pop_sel_pop, pos2 = pop_sel_pop)
    pairs_pop$type = ifelse(pairs_pop$pos1==pairs_pop$pos2, "within", "between")
    # pairs_pop$join_char = paste0(pairs_pop$pos1, "_", pairs_pop$pos2)
    # uniq_pairs = unique(pairs_pop$join_char)
    # more_colors = colorRampPalette(col)(length(col)^2)
    if(length(col) != 2){ 
      print("ERROR, incorrect color length")
      return()
    } else {
      pairs_pop$color = ifelse(pairs_pop$type == "within", col[1], col[2])
      return(list(pairs_pop$color, col))
    }
}


plotRegression = function(pop_sel, pdist, gdist, color, label, pch = 16, alpha = 1){
    # matrix to vector
    sel = as.vector(lower.tri(gdist))
    gdist = as.vector(gdist)[sel]
    pdist = as.vector(pdist)[sel]
    color = col2alpha(color, alpha)
    plot(gdist, pdist, pch = pch, col=color, cex = 0.5, xlab = "Genetic Relationship", ylab = "Phenotype distance", main = toString(pop_sel))
    if(length(pop_sel) > 1){
      legend("topright", legend = c("within", "between"), fill = label)
    }
    #reg = lm(pdist~gdist)
    #abline(reg)
}


col2alpha <- function(col, alpha){
  res = lapply(col, function(x){
    rgb(col2rgb(x)[1,1]/255, col2rgb(x)[2,1]/255, col2rgb(x)[3,1]/255, alpha)
  })
  res = unlist(res)
}

