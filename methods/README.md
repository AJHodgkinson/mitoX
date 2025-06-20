# 1. Alignment of RNA sequencing data:

For GTEx data: for each sample the BAM file was downloaded, converted to fastq, aligned to a reference genome with STAR (including generating gene count information), reads were filtered (as per methods), and modification and cleavage rate estimates were generated.  See alignment directory for scripts used for each step.

# 2. eQTL Generation:

An infosheet was generated for each tissue type, containing information on DNA sample name, RNA sample name and any covariates being used for eQTL mapping (various models used). mitogwas (https://github.com/AJHodgkinson/mitogwas) was then used for each tissue type to perform eQTL mapping with plink.

# 3. Generate Genetic Scores for mtDNA-encoded transcript abundance

Genetic, expression and covariate data was used for each tissue types that was generated in the previous step.  mitoimpute () was then used to format files correctly for FUSION, to run FUSION and then to test models against real data.
