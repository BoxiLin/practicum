#' Subseting individuals in phenotype and genotype data ####
#'
#' Subset the genetics and phenotype data by their "sample.ID", remove all missed phenotype values.
#'
#' @param input_rds Path to .rds of FBM genetic data
#' @param out_put_bed Path of output aligned genetic data
#' @param pheno Phenotype data, contain "sample.ID" column
#'
#' @return .plink output
qc_subset <- function(input_rds, output_bed, pheno) {
  psa <- snp_attach(input_rds)
  phenotype <- left_join(psa$fam, pheno, by="sample.ID")

  message("Subseting bigSNP for pa...", Sys.time())
  subset(psa, backingfile = bigsnpr::sub_bed(output_bed),
         ind.row = rows_along(psa$fam)[!is.na(phenotype$pa_age)])
  message("Writting .bed for pa...", Sys.time())
  x <- snp_attach(bigsnpr::sub_bed(output_bed, "rds"))
  snp_writeBed(x, bedfile = output_bed)
  message("Done!", Sys.time())
}

