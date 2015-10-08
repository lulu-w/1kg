
#############
## Process Results from
##
## $plink2 
##
#################


# Get a list of samples from All.chr.vcf.gz
bcftools query -l All.chr.vcf.gz > All.samples.txt

# Subset All.chr.vcf.gz by selected population samples
bcftools view -S asn.txt All.chr.vcf.gz -O z -o All.chr.asn.vcf.gz

# Calculate genetic relationship matrix using plink2
# All population
plink2 --vcf All.chr.vcf.gz --maf 0.05 --make-rel square gz --out all.chr.grm
# Eastern Asian populations only
plink2 --vcf All.chr.asn.vcf.gz --maf 0.05 --make-rel square gz --out all.chr.asn.grm

##################################
##
## R code
##
##################################

library(data.table)
options(stringsAsFactors = F)

####################
# Read data
####################

# Read genetic relationship matrix (GRM)
grm = fread("./bcf_concat/all.chr.asn.grm.rel")
grm_id = read.delim("bcf_concat/all.chr.asn.grm.rel.id", header = F)


# Read TE count as phenotype and calculate phenotype distance
teCount = read.delim("./abilbeast_input/presenceExp.txt", header = T)
sampleInfo = read.delim("./abilbeast_input/sampleInfo.txt", header = T)

#######################
# File pre-processing
#######################

# sample info and counts
sampleInfoList = split(sampleInfo, f = as.factor(sampleInfo$group))
subpopInfoList = split(sampleInfo, f = as.factor(sampleInfo$population))

tePopCountList = lapply(sampleInfoList, function(d){
	counts = unlist(teCount[1, d$id])
	counts = counts[order(names(counts))]
})

teSubpopCountList = lapply(subpopInfoList, function(d){
	counts = unlist(teCount[1, d$id])
	counts = counts[order(names(counts))]
})



popCounts = tePopCountList[["EAS"]]

subpopCounts = teSubpopCountList[["CHB"]]

# Genetic distances
grm_id = grm_id$V1
setnames(grm, grm_id)
grm = as.matrix(grm[, order(names(grm)), with=F])

subpopGrmList = lapply(subpopInfoList, function(d){
	colNum = which(colnames(grm) %in% d$id)
	subGrm = grm[colNum, colNum]
})

subpopGrm = subpopGrmList[["CHB"]]
sel = as.vector(!upper.tri(subpopGrm))
subpopGrm = as.vector(subpopGrm)
subpopGrm = subpopGrm[sel]

delta_y = as.matrix(dist(as.matrix(subpopCounts), diag = T))
delta_y = as.vector(delta_y)[sel]

plot(subpopGrm, delta_y, pch = ".", col = rgb(1, 0, 1, 0.5)); dev.off()



sel = as.vector(!upper.tri(grm))
grm = as.vector(grm)
grm = grm[sel]

# Calculate phenotype distance delta_y
delta_y = as.matrix(dist(as.matrix(popCounts), diag = T))
delta_y = as.vector(delta_y)[sel]




