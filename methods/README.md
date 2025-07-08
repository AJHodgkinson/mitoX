The following is an overview of the steps undertaken to analyse the role of mtDNA-encoded transcript abundance in common idsease risk. For more detailed codes and descriptions, use the subfolders above.

# 1. Alignment of GTEx RNA sequencing data:

For each sample the BAM file was downloaded, converted to fastq, aligned to a reference genome with STAR (including generating gene count information), reads were filtered (as per methods), and modification and cleavage rate estimates were generated.  See alignment directory for scripts used for each step.

# 2. eQTL Generation:

An infosheet was generated for each tissue type, containing information on DNA sample name, RNA sample name and any covariates being used for eQTL mapping (various models used). Mitogwas (https://github.com/AJHodgkinson/mitogwas) was then used for each tissue type to perform eQTL mapping with plink.

# 3. Generation of Genetic Scores for mtDNA-encoded transcript abundance

Genetic, expression and covariate data was used for each tissue type and gene that was generated in the previous step.  Mitoimpute () was then used to format files correctly for FUSION, to run FUSION and then to test models against real data.

# 4. Imputation of mtDNA-encoded transcript abundance into UKBB participants and tests against disease and quantitative traits


