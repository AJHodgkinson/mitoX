# Imputation

Using genetic score models generated using mitoimpute, mtDNA-encoded transcript abundances are imputed into each participant in UKBB for 15 genes across 49 tissue types, using different models filtered at different P-values (1,1e-1,1e-5) and using different machine learning prediction models (lasso, enet, blup):

```
cd /scratch/grp/hodgkinsonlab/shared/GTEx_Data/GTEx_Models
for i in *RDat; do stub=$(echo $i | awk -F "/" '{ print $1}');
perl impute_only_local_threads.pl ${stub}.enet ${stub}.wgt.RDat enet ukb.imp.maf1.miss5.hwe.MT.removed hg19tohg20 $stub/results/${stub}.selected.bim;
perl impute_only_local_threads.pl ${stub}.enet ${stub}.wgt.RDat lasso ukb.imp.maf1.miss5.hwe.MT.removed hg19tohg20 $stub/results/${stub}.selected.bim;
perl impute_only_local_threads.pl ${stub}.enet ${stub}.wgt.RDat blup ukb.imp.maf1.miss5.hwe.MT.removed hg19tohg20 $stub/results/${stub}.selected.bim;
done;
```
