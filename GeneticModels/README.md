# Genetic Models

Genetic score files are available for 15 genes encoded in mtDNA across 49 different tissue types.  For each gene/tissue combination, there are 9 different models to choose form (3 x P value thresholding for selected SNPs, 3 x ML method for score generation). 

To generate predicted transcript abundance for your samples (where you have genetic data available, use:

```plink2 --bfile <YOUR_GENETIC_DATA> --score <SCORE_FILE> --out <OUT_FILE_NAME>```

Where:

```<YOUR_GENETIC_DATA>``` is the genetic data for your samples in PLINK binary format

```<SCORE_FILE>``` is the score file model you want to use (gene, tissue, model)

```<OUT_FILE_NAME>``` is the name that you want to use for the output file

This approach averages scores across genetic loci per sample, to sum scores instead use:

```plink2 --bfile <YOUR_GENETIC_DATA> --score <SCORE_FILE> cols=+scoresums --out <OUT_FILE_NAME>```

Data is available at: https://figshare.com/projects/Nuclear_genetic_modulation_of_tissue-specific_mitochondrial_RNA_processing_contributes_to_common_disease_risk/256064
