# eQTL mapping

Infosheets were generated for each tisssue type (downloadable within this section), detailing DNA and RNA IDs, and any covarites to be used in eQTL mapping. Mitogwas (https://github.com/AJHodgkinson/mitogwas) was then used for each tissue type to perform eQTL mapping with plink, automatically creating PEER factors and genetic principle components to use as covariates, along with any additional covariates specified in the infosheet.  Mitogwas was the run for each tissue type seperately as follows:

```nextflow run main.nf --rnaDir "RNAseqDIR" --bed "BEDFILE" --bim "BIMFILE" --fam "FAMFILE" --infoSheet "INFOSHEET" --gtfFile "GTFFILE" --dataName "NAME" -profile singularity/docker```

This was run with both full and European only samples.


