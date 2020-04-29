library("bigsnpr")

bedfile <- "../backingfiles/canadian_pa.bed"
bedout <- "../backingfiles/canadian_pa_norel.bed"

# PATH to plink 2
plink2 <- "/hpf/tools/centos6/plink/2.0/plink2"

print(paste("# of cores", nb_cores()))

print("Start..............")

# Create new .bed file: "backingfiles/canadian_pa_norel.bed"

snp_plinkKINGQC(plink2,
                bedfile.in = bedfile,
                bedfile.out = bedout,
                make.bed = T,
                ncores = 12)

message("finish bedding!!")

df_rel <- snp_plinkKINGQC(plink2,
                          bedfile.in = bedfile,
                          bedfile.out = bedout,
                          make.bed = F,
                          ncores = 12)

saveRDS(df_rel, file = "../result/kinship_relation.rds")

# df_rel <- readRDS("../result/kinship_relation.rds")
# 
# hist(df_rel$KINSHIP, breaks = 100)
# 
# abline(v = c(0.0884, 0.125, 0.25, 0.5), col = "red", lty = c(1,2,2,2))
# 
# print("nice!!!!!!! no error")
