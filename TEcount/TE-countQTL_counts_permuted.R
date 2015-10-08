teCountQTL <- function(chunkName){
  
##########################
#
# Function to generate genotype permuted eQTL association results
#
##########################

options(stringsAsFactors = FALSE)
source("functions.R")
require(permute)

SNP_file_name_allPop = paste0("./abilbeast_input/", chunkName, ".genotype")
snps_location_file_name = paste0("./abilbeast_input/", chunkName, ".snploc")
expression_file_name = "./abilbeast_input/presenceExp.txt"
gene_location_file_name = "./abilbeast_input/presenceLoc.txt"
#covariates_file_name = "./abilbeast_input/covariates.txt"

snps <- read.table(SNP_file_name_allPop, header = T, stringsAsFactors = F)
snploc <- read.table(snps_location_file_name, header = T, stringsAsFactors = F)
snploc = snploc[!duplicated(snploc$snp), ]
geneExp <- read.table(expression_file_name, header = T, stringsAsFactors = F)
geneloc <- read.table(gene_location_file_name, header = T, stringsAsFactors = F)
#covar <- read.table(covariates_file_name, header = T, stringsAsFactors = F)
sample_info = read.table("./abilbeast_input/sampleInfo.txt", header = T)
colnames(snps)[1] = tolower(colnames(snps)[1])


######################
# Covariates

sample_info_pop = split(sample_info, f=sample_info$group)
lapply(sample_info_pop, function(d){
    xgroup = d$group[1]
    x = d[,2:3]
    rownames(x) = d$id
    x$population = as.numeric(as.factor(x$population))
    x$gender = as.numeric(as.factor(x$gender))
    x = t(x)
    x = data.frame(id=rownames(x), x)
    write.table(x, paste0("./temp/covariates_", xgroup, "_", chunkName, ".txt"), quote = F, sep = "\t", row.names = F)
})
pop_individuals = lapply(sample_info_pop, function(d){ d$id })


#####################
# SNP genotype

pop_snps = lapply(pop_individuals, function(x)snps[, c("id", x)] )
af_cut = 0.05
pop_snps_filtered = lapply(pop_snps, function(d){
    af_lb = af_cut * (ncol(d)-1) * 2
    x = d[which(rowSums(d[, -1]) > af_lb), ]
})
mapply(write.table, pop_snps_filtered, 
       file = paste0("./temp/snps_", names(pop_snps_filtered), "_", chunkName, ".txt"), 
       quote = rep(FALSE, length(pop_snps_filtered)), 
       row.names = rep(FALSE, length(pop_snps_filtered)), 
       sep = rep("\t", length(pop_snps_filtered)) )

# SNP number log

chunk_snp_count=unlist(lapply(pop_snps_filtered, nrow))
chunk_snp_count=c(chunkName, chunk_snp_count)
write.table(t(as.matrix(chunk_snp_count)), "./abilbeast_output/chunk_snp_count.txt", quote=F, col.names=F, row.names=F, append=T)


#######################
# SNP locations

filtered_snps = lapply(pop_snps_filtered, function(x) x$id )
pop_snploc = lapply(filtered_snps, function(x){ df = snploc
                                            rownames(df) = snploc$snp
                                            df[x, ]})
mapply(write.table, pop_snploc, 
       file = paste0("./temp/snploc_", names(pop_snploc), "_", chunkName, ".txt"), 
       quote = rep(FALSE, length(pop_snploc)), 
       row.names = rep(FALSE, length(pop_snploc)), 
       sep = rep("\t", length(pop_snploc)) )

###################
# Expression

pop_geneExp = lapply(pop_individuals, function(x) geneExp[, c("id", x)] )

# Shuffle
set.seed(12345) 
pop_geneExp_shuffled = lapply(pop_geneExp, function(d){
  ori_colnames=colnames(d)
  toShuffle=2:ncol(d)
  shuffled=toShuffle[shuffle(length(toShuffle))]
  newD=d[, c(1,shuffled)]
  colnames(newD)=ori_colnames
  return(newD)
})


mapply(write.table, pop_geneExp_shuffled, 
       file = paste0("./temp/geneExp_shuffled_", names(pop_geneExp), "_", chunkName, ".txt"), 
       quote = rep(FALSE, length(pop_geneExp)), 
       row.names = rep(FALSE, length(pop_geneExp)), 
       sep = rep("\t", length(pop_geneExp)) )

################################
# population Specific Results

# Global par

p.cut = 1e-2
popme = list()

for(i in 1:length(pop_geneExp)){
    xgroup = names(pop_geneExp)[i]
    covariates_file_name = paste0("./temp/covariates_", xgroup, "_", chunkName, ".txt")
    SNP_file_name = paste0("./temp/snps_", xgroup, "_", chunkName, ".txt")
    snps_location_file_name = paste0("./temp/snploc_", xgroup, "_", chunkName, ".txt")
    expression_file_name = paste0("./temp/geneExp_shuffled_", xgroup, "_", chunkName, ".txt")
    xsnps = read.table(snps_location_file_name)
    if(nrow(xsnps)==1){
        popme[[i]] = NULL
    }else{
      res = runME( SNP_file_name, 
                   snps_location_file_name, 
                   expression_file_name, 
                   gene_location_file_name, 
                   covariates_file_name, 
                   output_file_name = "./all_out.txt", 
                   useModel = modelLINEAR, 
                   par.plot = 'qqplot', 
                   verbose = F, 
                   pvOutputThreshold = p.cut, 
                   noFDRsaveMemory = F)
      popme[[i]] = res$all$eqtls
    }
    file.remove(covariates_file_name)
    file.remove(SNP_file_name)
    file.remove(snps_location_file_name)
    file.remove(expression_file_name)
}



#######################
# Process matrixEQTL results


if(length(popme)>0){
  
  # Manhattan plot pvalues
  
    require(dplyr)
    popme_pval = lapply(popme, function(d){
      if(is.null(d)){
        x = NULL
      }else{
        # d = d[order(d$snps, d$pvalue), c("snps", "pvalue")]
        d = d[order(d$snps, d$pvalue), ]
        x = d[!duplicated(d$snps), ]
      }
    })
    
    for(i in 1:length(popme_pval)){
      colnames(popme_pval[[i]])[1] = "snp"
      popme_pval[[i]] = left_join(popme_pval[[i]], pop_snploc[[i]], by = "snp")
      popme_pval[[i]] = popme_pval[[i]][order(popme_pval[[i]]$pos), ]
    }
    
    mapply(write.table, popme_pval, 
           file = paste0("./abilbeast_output/", chunkName, "_", names(pop_geneExp), ".txt"), 
           quote = rep(FALSE, length(popme_pval)), 
           row.names = rep(FALSE, length(popme_pval)), 
           sep = rep("\t", length(popme_pval)) )
    
    # Q-Q plot pvalues
    mapply(write.table, popme, 
           file = paste0("./abilbeast_output/", chunkName, "_AllPvalues_", names(pop_geneExp), ".txt"), 
           quote = rep(FALSE, length(popme_pval)), 
           row.names = rep(FALSE, length(popme_pval)), 
           sep = rep("\t", length(popme_pval)) )
    
    print(paste("Chunk", chunkName, "done!"))
    
    return()
}else{
    print(paste("Chunk", chunkName, "has no significant association!"))  
    return()
}

}
