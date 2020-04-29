#' Clean phenotype dataset
#'  1. assign CFTR score
#'  2. modify variables for data and structures visualization
#'  3. filter invalid observation, extract response and covariates for gwas; 
#'  
#' Output: "phenotype_gwas" outcome for gwas analysis
#'    
#'
#'  
#'  #######################################################

library(dplyr)
#library(readxl)
#library(lubridate)
#library(eeptools)
#library(GGally)
#library(boxi)

### read data ####

# Filter out US samples
raw1 = readr::read_csv("../raw/CFGeneModifierStudy-Boxirace_DATA_2020-04-03_1458.csv") %>% dplyr::filter(substr(individual_name,1,1) != "X")
raw2 = readr::read_csv("../raw/Psa_2752genotyped_Mar2_2020.csv")
raw <- full_join(raw1, raw2, by = "individual_name")

dplyr::glimpse(raw)

### phenotype data frame  ######
phenotype <- raw %>% 
  filter(!(is.na(pac_age.x) & is.na(pac_age.y))) %>%
  filter(panc_function==0) %>%
  mutate(pa = factor(pa.x,levels = c(-999, 0:3), 
                labels = c("Missing", "No growth", 
                           "Grew once", "Sporadic growth",
                           "Chronically infected")),
         sex = as.factor(sex),
         ethnic = factor(ethnic_orig,levels = 0:4, 
                         labels = c("european", "asian", 
                                    "african", "other",
                                    "mixed")),
         cftr = boxi::cftr(allele1, allele2))
phenotype$pac_age.y[is.na(phenotype$pac_age.y)] <- phenotype$pac_age.x[is.na(phenotype$pac_age.y)]
phenotype <- phenotype %>%
  select(individual_name, site,pac_date.y, pac_age.y, pa_source, dob, sex, pa, ethnic, cftr)

colnames(phenotype)[colnames(phenotype) == "individual_name"] <- "sample.ID"
colnames(phenotype)[colnames(phenotype) == "pac_date.y"] <- "pac_date"
colnames(phenotype)[colnames(phenotype) == "pac_age.y"] <- "pac_age"

### update sex ####
phenotype[phenotype$sample.ID=="BSP0080-03","sex"] <- "F"
phenotype[phenotype$sample.ID=="BSP0049-03","sex"] <- "M"
phenotype[phenotype$sample.ID=="AFO3932-03","sex"] <- "M"
phenotype[phenotype$sample.ID=="ACH1202-03","sex"] <- "M"

sum(is.na(phenotype$pac_age))


saveRDS(phenotype, "../clean/phenotype_pac.rds")
