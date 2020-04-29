library(dplyr)

### 1. preprocess ####
gwas <- readRDS("output/assoc_pa.rds")

# the following chunk run once
# epical <- data.table::fread("4016snps_3988apical_38_ATP12A_PRSS1_SLC6A14.txt")[,c(1,2,4:9)]
# colnames(epical)[1:3] <- c("Chr", "variant.id", "Bp")
# saveRDS(epical,"temp/epical.rds")


irf1 <- readRDS("../post/irf1.rds")
irf1_set <- dplyr::inner_join(irf1, gwas, by = "variant.id") %>% filter(!is.na(Score.pval))
hist(irf1_set$Score.pval, breaks = 100)

epical <- readRDS("temp/epical.rds")
epical_set <- dplyr::inner_join(epical, gwas, by = "variant.id") %>% filter(freq > 0.05 & hwe_p > 1e-10)
hist(epical_set$Score.pval, breaks = 100)


epical_snp <- snpStats::read.plink(bed = "plink.bed",
                         bim = "plink.bim",
                         fam = "plink.fam")
cor(as.matrix(epical_snp$genotypes[1:10, 1:10]))

### 2. Cauchy combination test ####
cauchy_test <- function(p, w = 1/length(p)) {
  T = sum(w*tan((0.5-p)*pi))
  p = 2*(1-pcauchy(abs(T)))
  return(c("T"=T, "p_value"=p))
}
cauchy_test(epical_set$Score.pval)
cauchy_test(irf1_set$Score.pval)


### 3. Set-based test ####

set_based_test<-function(summary_stats,LDmatrix,total_number_genes,alpha){
  statistic<-sum(summary_stats^2)
  m=length(summary_stats)
  matrix_mid<-LDmatrix
  eig_result<-eigen(matrix_mid)
  eigenvalues<-eig_result$values
  pv<-abs(imhof(statistic,eigenvalues,h=rep(1,m),delta=rep(0,m))$Qq)
  ##using bonferoni at the first stage
  if (pv<(alpha/total_number_genes)){
    result=TRUE
  }else{
    result=FALSE
  }
  return(result)
}



write.table(epical_set$variant.id, file ="epical_snp_id.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)



ld <- data.table::fread("sampled_ld.ld")
length(unique(ld$SNP_A))
length(unique(ld$SNP_B))
length(unique(c(ld$SNP_B, ld$SNP_A)))


matrix <- diag(ncol =3565, nrow = 3565)
colnames(matrix) = epical_set$variant.id
rownames(matrix) = epical_set$variant.id


for (i in 1:9210) {
  matrix[ld$SNP_A[i],ld$SNP_B[i]] <- ld$R2[i]
}



for (i in 1:3565) {
  for (j in i:3565) {
    matrix[j,i] = matrix[i,j]
  }
}


set_based_test(summary_stats = epical_set$Score.Stat,
               LDmatrix = ldm, total_number_genes = 1,alpha = 0.05)


