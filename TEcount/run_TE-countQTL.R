
chunkList = read.table("chunkList_chr1.txt", header = T)
chunkList = paste0(chunkList$chr,"_",chunkList$starts,"_",chunkList$ends)


#chunkList = c("1_1_1000000", "1_1000001_2000000")


source("TE-countQTL.R")

start_time = proc.time()

lapply(chunkList, teCountQTL)

elapsed = proc.time() - start_time
print(elapsed)
