library(bigsnpr)
library(dplyr)

message("making bk file...")
psa <- bigsnpr::snp_readBed2(
  bedfile = "../backingfiles/canadian_pa_norel.bed",
  backingfile = "../backingfiles/canadian_pa_norel")
message("Done making bk file...")

message("Read RDS: contain no missing Pa...", Sys.time())
psa <- snp_attach("../backingfiles/canadian_pa_norel.rds")

phenotype = readRDS("../../phenotype/clean/phenotype_eda.Rds")
phenotype <- left_join(psa$fam, phenotype, by="sample.ID")
saveRDS(phenotype, file = "../backingfiles/phenotype_pa.RDS")

y = phenotype$pa_age
message(paste("missed phenotype:",sum(is.na(y))))

sex <- 2-as.matrix((as.numeric(phenotype$sex)))
dob <- scale(as.numeric(phenotype$dob))
cftr <- phenotype$cftr
pca <- readRDS("../result/pca_april3.RDS")

cov <- as.matrix(sex)

message("Run gwas......", Sys.time())
gwas <- big_univLinReg(psa$genotypes,
                       y.train = y,
#                       covar.train = cov,
                       ncores = 10)
message("Done gwas.....", Sys.time())
saveRDS(gwas, file = paste("../result/gwas/gwas_",Sys.Date(),".RDS", sep = ""))

message("GWAS Success!!!!!!!!!!!!!!!!!!!!!!!!!!", Sys.time())

