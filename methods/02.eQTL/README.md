mitogwas (https://github.com/AJHodgkinson/mitogwas) was then used for each tissue type to perform eQTL mapping with plink:

nextflow run main.nf --rnaDir "RNAseqDIR" --bed "BEDFILE" --bim "BIMFILE" --fam "FAMFILE" --infoSheet "INFOSHEET" --gtfFile "GTFFILE" --dataName "NAME" -profile singularity/docker

This was run with both full and European only samples.


