library(bigsnpr)
library(dplyr)


tmpfile <- "backingfiles/canadian_merged_norel_QC"
source_bed <- "../../../Data_query/vcf2bed/non_fam_id/canadian_merged_norel"

plink <- "/hpf/tools/centos6/plink/1.90b3x/plink"
qc <- snp_plinkQC(plink.path = plink,
                    prefix.in = source_bed,
                    prefix.out = tmpfile,
                    file.type = "--bfile",  # the default (for ".bed")
                    maf = 0.01,
                    geno = 0.05,
                    mind = 0.05)


rdsfile <- snp_readBed2(qc, backingfile = tmpfile)