library(SNPRelate)
library(GWASTools)
library(GENESIS)
library(dplyr)
library(hexbin)

# # 1. PC-AiR
# pheno <- readRDS("phenotype.rds")
# 
# mypcair <- readRDS("temp/mypcair.rds")
# pallete <- c("red", "purple", "blue", "grey", "green")
# names(pallete)<- levels(pheno$ethnic)
# pallete
# 
# par(mfrow = c(1,3))
# cols <- pallete[pheno$ethnic]
# for (i in c(1,3,5)) {
#   plot(mypcair,vx = i+1, vy = i, col = scales::alpha(cols,0.7))
#   legend(x = "topleft", legend = names(pallete), fill = pallete)
# }
# 
#
#
# # 2. PC-Relate
# par(mfrow = c(1,1))
# plot(mypcair$values, type = "o", ylab = "Singular value", main = c("Scree plot"))
# 
#      
# mypcrelate <- readRDS("temp/mypcrelate.rds")
# plot(mypcrelate$kinBtwn$k0, mypcrelate$kinBtwn$kin, xlab="k0", ylab="kinship")
#
#
#
# 5. Post GWAS analysis
#
# Following only run once
# gwas <- readRDS("output/assoc.rds")
# hwe <- readRDS("../qc/hwe.rds")
# colnames(hwe)[colnames(hwe)=="SNP"] <- "variant.id"
# gwas <- dplyr::left_join(gwas, hwe, by = "variant.id")
# gwas_filtered = gwas[,c("variant.id", "chr","pos", "n.obs" , 
#                         "freq", "Score", "Score.SE", "Score.Stat", 
#                         "Score.pval", "A1" , "A2", "GENO", "P")]
# colnames(gwas_filtered)[colnames(gwas_filtered)=="P"]<- "hwe_p"
# saveRDS(gwas_filtered, file = "output/assoc_pa.rds")


nullmod <-readRDS("output/nullmod.rds")
varCompCI(nullmod, prop = TRUE)

assoc <-readRDS("output/assoc_pa.rds")
assoc$hwe_p[assoc$chr=="X"]=1

assoc_plot <- dplyr::filter(assoc, Score.pval <0.01 & freq>0.05 & hwe_p > 1e-10 & !is.na(Score.pval))[1:9]

colnames(assoc_plot)<- c("SNP","CHR","BP", "N", "freq", "zscore","score.se", "score.st", "P")
assoc_plot$CHR[assoc_plot$CHR=="X"] = "23"
assoc_plot$CHR <- as.numeric(assoc_plot$CHR)
qqman::manhattan(assoc_plot,ylim = c(2,8), annotatePval = 0.00001,
                 main = "Genome-wide Manhattan plot of associations with Pa age")

log10p_value <- data.frame(lmm = -log10(assoc$Score.pval), 
                           lm = -log10(2*(1-pnorm(abs(assoc2$score))))) %>% 
  filter(!is.na(lmm) & !is.na(lm) & assoc$freq>0.05 & assoc$P > 1e-10 )

hexbinplot(lm ~ lmm, data = log10p_value, aspect = '1', xbins = 300,
           xlim = c(0, 7), ylim = c(0, 7),xlab = "-log10(p) from Linear mixed model",
           ylab =  "-log10(p) from Linear model", 
           panel = function(x, y, ...) {
             panel.hexbinplot(x, y, ...)
             lattice::panel.abline(a = 0, b = 1)})
