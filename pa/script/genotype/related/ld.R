library(snpStats)
# id <- readRDS("irf1_id.rds")
# irf1 <- read.plink(bed = "../backingfiles/canadian_pac_all.bed",
#                    bim = "../backingfiles/canadian_pac_all.bim",
#                    fam = "../backingfiles/canadian_pac_all.fam")
# 
# irf1set <- irf1$genotypes[,id]
# saveRDS(irf1set, file = "irf1_matrix.rds")
# 
# message(paste("Calculating LD...", Sys.time()))
# ldm <- ld(x=irf1set, depth = 500, stats = "Covar")
# 
# saveRDS(ldm, file = "ld_matrix.rds")
# 
# message(paste("Done LD!!!", Sys.time()))


ld <- readRDS("ld_matrix.rds")
pdf(file = "plot.pdf")
image(ld)
dev.off()
message("done")
# for(i in 1:77295){
#   ld[i,1]<- 1
#   print(i)
# }
