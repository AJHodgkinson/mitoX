## Genetic Models

# Model Descriptions

Genetic score files are available for 15 genes encoded in mtDNA across 49 different tissue types.

Data is available at: https://figshare.com/projects/Nuclear_genetic_modulation_of_tissue-specific_mitochondrial_RNA_processing_contributes_to_common_disease_risk/256064

Each score file contains the following columns:

1. SNP ID (chr_position, in genome build stated in the file name)
2. Effect allele (the allele used for scoring)
3. Weight (effect size used in the genetic score)

For each gene-tissue pair, 9 models are provided:
- 3 SNP inclusion thresholds (based on association P-values)
- 3 modelling approaches (e.g. elastic net, etc.)

Models where mtDNA-encoded transcript abundance predictions are associated with traits are shown in Supplementary Tables 8 and 10.

Note: It is important to ensure that alleles in your dataset are aligned to the effect allele in the score file. Mismatches or strand flips may lead to incorrect predictions.

# Model Usage

To generate predicted transcript abundance for your samples (where you have genetic data available, use:

```plink2 --bfile <YOUR_GENETIC_DATA> --score <SCORE_FILE> --out <OUT_FILE_NAME>```

Where:

```<YOUR_GENETIC_DATA>``` is the genetic data for your samples in PLINK binary format

```<SCORE_FILE>``` is the score file model you want to use (gene, tissue, model)

```<OUT_FILE_NAME>``` is the name that you want to use for the output file

This approach averages scores across genetic loci per sample, to sum scores instead use:

```plink2 --bfile <YOUR_GENETIC_DATA> --score <SCORE_FILE> cols=+scoresums --out <OUT_FILE_NAME>```

# Model Example

Using the files available in this repository, running the following command should give idenitical output as found in "test.out.sscores":

```plink2 --bfile 1000G_chr22_50inds --score WholeBlood.ENSG00000228253.1.blup.1e1.weights.hg19.txt cols=+scoresums --out test.out```




