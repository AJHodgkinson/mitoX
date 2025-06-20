The methods described here are to accompany the study "Nuclear-genetic modulation of tissue-specific mitochondrial RNA processing contributes to common disease risk", published as...

# 1. Alignment of RNA sequencing data:

For GTEx data: for each sample the BAM file was downloaded, converted to fastq, aligned to a reference genome with STAR (including generating gene count information), reads were filtered (as per methods), and modification and cleavage rate estimates were generated.  See alignment directory for scripts used for each step.

# 2. eQTL Generation:

An infosheet was generated for each tissue type, containing information on DNA sample name, RNA sample name and any covariates being used for eQTL mapping (various models used). mitogwas (https://github.com/AJHodgkinson/mitogwas) was then used for each tissue type to perform eQTL mapping with plink.
