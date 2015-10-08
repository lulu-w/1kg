#chunkList = c("1_1_1000000", "1_1000001_2000000")

chunkList = read.table("chunkList_chr1.txt", header = T)
chunkList = paste0(chunkList$chr, "_", chunkList$starts, "_", chunkList$ends)

library(parallel)

start_time = proc.time()


source("TE-countQTL.R")
mclapply(chunkList, teCountQTL,  mc.cores = 10)


elapsed = proc.time() - start_time
print(elapsed)