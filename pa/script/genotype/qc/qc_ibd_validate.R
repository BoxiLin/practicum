library("bigsnpr")
library("dplyr")


bedfile <- "../backingfiles/canadian_pa.bed"

# PATH to plink 
plink <- "/hpf/tools/centos6/plink/1.90b3x/plink"

message("IBD Start..............", Sys.time())

# check relationship pairs
df_rel <- snp_plinkIBDQC(plink,
                         bedfile.in = bedfile,
                         do.blind.QC = FALSE,
                         ncores = 12)

saveRDS(df_rel, file = "../result/ibd_relation.rds")

relationship <- function(z0, z1){
  z2 = 1- z0-z1
  rel<- ifelse(z0>0.3 & z2<0.1, "UNRELATED",
               ifelse(z0>0.1 & z1>0.3 & z2>0.1, "SIB-PAIRS",
                      ifelse(z0<0.1 & z1>0.8&z2<0.2,"PARENT-OFFSPRING",
                             ifelse(z2>0.9,"TWINS/DUPLICATES","UNDEFINED"))))
}

### Analysing Relatedness
library(dplyr)
df_rel <- readRDS("../result/ibd_relation.rds")
tf_rel <- data.table::fread("psa_truffle.ibd", header = TRUE)
tf_rel$pi <- 0.5*tf_rel$IBD1 + tf_rel$IBD2
tf_pair <- tf_rel %>% dplyr::filter(pi >= 0.08)
tf_pair$rel <- relationship(tf_pair$IBD0, tf_pair$IBD1)
df_rel$rel <- relationship(df_rel$Z0,df_rel$Z1)


library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(as.factor(df_rel$rel))
colScale <- scale_colour_manual(name = "Relatedness",values = myColors)

library(ggplot2)
p1 <-qplot(IBD0, IBD1, data = tf_pair,xlim = c(0,1),
           ylim = c(0,1), 
      color=rel, alpha=0.5)+colScale+
  ggtitle("TRUFFLE estimated IBD, 196 pairs, pi_hat>0.08")


p2<-qplot(Z0, Z1, data = df_rel,xlim = c(0,1),
          ylim = c(0,1), 
      color=rel, alpha=0.5)+colScale+
  ggtitle("PLINK estimated IBD, 182 pairs, pi_hat>0.08")


gridExtra::grid.arrange(p1,p2, ncol = 2)
