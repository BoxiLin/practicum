#' This script conducted association test for Pa age with 2740 Canadian sample,  based on GENESIS 
#' (https://www.bioconductor.org/packages/release/bioc/html/GENESIS.html) package.
#' 
#' - Reference: 
#'     1. http://faculty.washington.edu/tathornt/SISG2015/lectures/Taipei2015/Taipei2015session06.pdf
#'     2. https://www.bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test.html
#' 
#' 
#' - input: 
#'     1. binary plink format (.bed, .bim, .fam) of genotype data after QC
#'     2. sample relatedness estimated from KING
#'     3. phenotype data
#'     
#' - output:
#'     1. outcome of null models
#'     2. summary statistics of association test
#'     

library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(dplyr)
library(hexbin)

### 1. Creating GDS file #### 
message(paste("Start convert to gds......", Sys.time()))
SNPRelate::snpgdsBED2GDS(bed.fn = "../backingfiles/canadian_pac_all.bed",
                         bim.fn = "../backingfiles/canadian_pac_all.bim",
                         fam.fn = "../backingfiles/canadian_pac_all.fam",
                         out.gdsfn = "../backingfiles/genotype_pac.gds")
message(paste("Done gds!!!!!!", Sys.time()))

### 2. Prunning for PCA ####
gds <- SNPRelate::snpgdsOpen( "../backingfiles/genotype_pa.gds")
message(paste("Start pruning......", Sys.time()))
snpset <- SNPRelate::snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
                                     ld.threshold=sqrt(0.1), verbose=TRUE,
                                     num.thread = 12)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
head(pruned)
saveRDS(pruned, file = "temp/pruned.rds")
snpgdsClose(gds)
message(paste("Done pruning......", Sys.time()))


### 3. Running PC-Air ####
message(paste("Start PC-Air......", Sys.time()))
# read in KINGship matrix
KINGmat <- kingToMatrix(c("king.kin0","king.kin"))
pruned <- readRDS("temp/pruned.rds")

psa_geno <- GdsGenotypeReader(filename = "../backingfiles/genotype_pac.gds")
psa_genoData <- GenotypeData(psa_geno)
psa_genoData

mypcair <- pcair(psa_genoData, kinobj = KINGmat, divobj = KINGmat,
                 snp.include = pruned)

saveRDS(mypcair, "temp/mypcair_pac.rds")

### 4. Running PC-Relate (optional) ####
message(paste("Start PC-Relate......", Sys.time()))
psa_genoData <- GenotypeBlockIterator(psa_genoData, snpInclude=pruned)
mypcrelate <- pcrelate(psa_genoData, pcs = mypcair$vectors[,1:5],
                       training.set = mypcair$unrels)
saveRDS(mypcrelate, "temp/mypcrelate_pac.rds")
message(paste("Done everything!!!!!!", Sys.time()))


### 5. Preprocessing for GWAS ####
rm(list=ls())

# output from step 3 & 4
mypcair <- readRDS("temp/mypcair_pac.rds")
mypcrel <- readRDS("temp/mypcrelate_pac.rds")

# read in phenotype and create phenotype matrix
phenotype <- readRDS("phenotype_pac.rds")
mydat <- data.frame(scanID = mypcair$sample.id,
                    pheno = phenotype$pac_age,
                    mypcair$vectors[,1:6],
                    dob = scale(as.numeric(phenotype$dob)),
                    sex = phenotype$sex.y,
                    cftr = phenotype$cftr)

scanAnnot <- ScanAnnotationDataFrame(mydat)
scanAnnot

psa_geno <- GdsGenotypeReader(filename = "../backingfiles/genotype_pac.gds")
psa_genoData <- GenotypeData(psa_geno, scanAnnot = scanAnnot)
psa_genoData

myGRM <- pcrelateToMatrix(mypcrel)

### 6. Fitting null model ####
nullmod <- fitNullModel(scanAnnot, outcome = "pheno",
                        covars = colnames(mydat)[3:11],
                        cov.mat = myGRM, family = gaussian)

### 7. SNP association testing ####
genoIterator <- GenotypeBlockIterator(psa_genoData, snpBlock = 10000)
assoc <- assocTestSingle(genoIterator, null.model = nullmod)

saveRDS(nullmod, file = "nullmod_pac.rds")
saveRDS(assoc, file = "assoc_pac.rds")

message(paste("Done everything!!!!!!", Sys.time()))