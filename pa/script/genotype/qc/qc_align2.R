library(bigsnpr)
library(dplyr)

### Align individuals in phenotype and genotype data ####

# Need to run only once
# getwd()
# 
# source_bed <- "../backingfiles/canadian_pa_norel.bed"
# tmpfile <- "../backingfiles/canadian_pa_norel"
# message("making bk file...")
# psa <- snp_readBed2(
#     bedfile = source_bed,
#     backingfile = tmpfile
# )
# message("done making bk file...")

psa <- snp_attach("../backingfiles/canadian_pa_norel.rds")
phenotype = readRDS("../../phenotype/clean/phenotype_eda.Rds")
phenotype <- left_join(psa$fam, phenotype, by="sample.ID")

saveRDS(phenotype,"../backingfiles/phenotype_pa.rds")

message("Done!")