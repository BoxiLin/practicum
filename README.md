# Appendix 6.3B

This is the online appendix for the code to reproduce my practicum project in Winter 2020, including programming code for phenotype data exploration, sample & genotype quality control, post-analysis plots and helper function package, etc.

The structure of the files is as the follows:
```{bash}
├── pkg
│   └── boxi
│       ├── DESCRIPTION
│       ├── NAMESPACE
│       ├── R
│       │   ├── cftr_score.R
│       │   ├── qc_subset.R
│       │   ├── qqunif_plot.R
│       │   └── txt2rds.R
│       ├── boxi
│       ├── boxi.Rproj
│       └── man
│           ├── qqunif.plot.Rd
│           └── txt2rds.Rd
├── reference
│   ├── gwas_implement
│   │   └── bigsnp.pdf
│   ├── pcombine
│   │   ├── cauchy.comb.test.pdf
│   │   └── combine\ p-values.pdf
│   └── phenotype
│       ├── cftr.pptx
│       ├── gwas.pdf
│       └── survival.pdf
└── script
    ├── genotype
    │   ├── backingfiles
    │   ├── gwas
    │   │   └── gwas.R
    │   ├── heritability
    │   │   └── submit.pbs
    │   ├── qc
    │   │   ├── qc_align2.R
    │   │   ├── qc_ibd_validate.R
    │   │   ├── qc_kinship.R
    │   │   ├── qc_pca_bed.R
    │   │   ├── qc_qc.R
    │   │   └── qc_subset.R
    │   ├── related
    │   │   ├── GENESIS.R
    │   │   ├── ld.R
    │   │   ├── output
    │   │   ├── post_plot.R
    │   │   ├── set_test.R
    │   │   └── temp
    │   └── submit.pbs
    └── phenotype
        ├── clean
        │   └── phenotype.R
        ├── eda
        │   ├── eda_phenotype.Rmd
        │   └── eda_phenotype.html
        ├── manipulate
        └── raw
```
