library(ggplot2)
library(dplyr)
library(bigsnpr)

### Run PCA #####
message("bedding........", Sys.time())
(obj.bed <- bed("../backingfiles/canadian_pa_all_ethics.bed"))

NCORES <- 11

message("bed PCAing.........", Sys.time())
obj.svd <- bigsnpr::bed_autoSVD(obj.bed, k = 20, 
                                ncores = NCORES)
message("saving.......", Sys.time())
saveRDS(obj.svd, file = "../result/pca_all_ethics.RDS")
message("done!  ", Sys.time())


## PCA analysis ####
obj.svd <- readRDS("../result/pca_all_ethics.RDS")
psa <- snp_attach("../backingfiles/canadian.rds")
phenotype = readRDS("../../phenotype/clean/phenotype_all_ethics.Rds")
phenotype <- left_join(psa$fam, phenotype, by="sample.ID") %>% filter(!is.na(pa_age))


plot(obj.svd)
p1 <- plot(obj.svd, type = "scores",scores = 1:2, coeff = 0.4) +
  aes(color = phenotype$ethnic, alpha = 0.4)
p2 <- plot(obj.svd, type = "scores",scores = 3:4, coeff = 0.4) +
  aes(color = phenotype$ethnic, alpha = 0.4)
gridExtra::grid.arrange(p1, p2, ncol=2)

plot(obj.svd, type = "loadings", loadings = 1:9, coeff = 0.4)

prob <- bigutilsr::prob_dist(obj.svd$u, ncores = 2)
S <- prob$dist.self / sqrt(prob$dist.nn)

ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)")

plot_grid(plotlist = lapply(1:6, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S) +
    scale_colour_viridis_c()
}), scale = 0.95)



