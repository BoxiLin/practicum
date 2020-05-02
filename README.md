# Appendix 6.3B

This is the online appendix for the code to reproduce my practicum project in Winter 2020, including programming code for phenotype data exploration, sample & genotype quality control, post-analysis plots and helper function package, etc.


## File structures

The structure of the files is as the follows:

```{bash}
├── pkg: Package of helper functions for GWAS, including data cleaning and visualization.
│   └── boxi
│       ├── DESCRIPTION
│       ├── NAMESPACE
│       ├── R
│       │   ├── cftr_score.R
│       │   ├── qc_subset.R
│       │   ├── qqunif_plot.R
│       │   └── txt2rds.R
│       ├── boxi
│       ├── boxi.Rproj
│       └── man
│           ├── qqunif.plot.Rd
│           └── txt2rds.Rd
├── reference
│   ├── gwas_implement
│   │   └── bigsnp.pdf
│   ├── pcombine
│   │   ├── cauchy.comb.test.pdf
│   │   └── combine\ p-values.pdf
│   └── phenotype
│       ├── cftr.pptx
│       ├── gwas.pdf
│       └── survival.pdf
└── script: main folder for GWAS
    ├── genotype: Dealing with genotypic data
    │   ├── backingfiles
    │   ├── gwas: baseline model
    │   │   └── gwas.R
    │   ├── heritability: heritability estimation
    │   │   └── submit.pbs
    │   ├── qc: genotypic quality control
    │   │   ├── qc_align2.R
    │   │   ├── qc_ibd_validate.R
    │   │   ├── qc_kinship.R
    │   │   ├── qc_pca_bed.R
    │   │   ├── qc_qc.R
    │   │   └── qc_subset.R
    │   ├── related: LMM model
    │   │   ├── GENESIS.R
    │   │   ├── ld.R
    │   │   ├── output
    │   │   ├── post_plot.R
    │   │   ├── set_test.R
    │   │   └── temp
    │   └── submit.pbs
    └── phenotype: phenotype data cleaning
        ├── clean
        │   └── phenotype.R
        ├── eda
        │   ├── eda_phenotype.Rmd
        │   └── eda_phenotype.html
        ├── manipulate
        └── raw
```



## QC steps with intermidiate input and output genotype data

0. 
   - raw genotypic data: 
     - `canadian_merged.bed/bim/fam`; 
   - script: 
     - `related/qc.R`




1. extract overlapped individuals in both genotype and phenotype data;

   - input: 
     - `canadian_merged.bed/bim/fam`
     -  `phenotype_pa.rds`
   - output: 
     - `canadian_pa_all.bed/bim/fam`

   


2. remove duplicates and twins
   - input:
     - `canadian_pa_all.bed/bim/fam`
   - output:
     - `canadian_pa_all_nodup.bed/bim/fam`
