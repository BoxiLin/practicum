library(bigsnpr)
library(dplyr)

### Align individuals in phenotype and genotype data ####

# Need to run only once
# getwd()
# 
# source_bed <- "../../../../Data_query/vcf2bed/non_fam_id/canadian_merged.bed"
# tmpfile <- "../backingfiles/canadian"
# message("making bk file...")
#  psa <- snp_readBed2(
#     bedfile = source_bed,
#     backingfile = tmpfile
# )
# message("done making bk file...")

psa <- snp_attach("../backingfiles/canadian.rds")
phenotype = readRDS("../../phenotype/clean/phenotype_pac.rds")
phenotype <- left_join(psa$fam, phenotype, by="sample.ID")

message("subseting bigSNP for pac...")
subset(psa, 
       ind.row = rows_along(psa$fam)[!is.na(phenotype$pac_age)],
       backingfile = "../backingfiles/canadian_pac_all")

message("writting .bed for pac...")
x <- snp_attach("../backingfiles/canadian_pac_all.rds")
snp_writeBed(x, bedfile = "../backingfiles/canadian_pac_all.bed")

message("Done!")


x <- snp_attach("../backingfiles/canadian_pac_all.rds")
phenotype <- left_join(x$fam, phenotype, by="sample.ID")
saveRDS(phenotype, "../related/phenotype_pac.rds")
