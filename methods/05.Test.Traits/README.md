# Test imputed mtDNA abundance against quantitative traits and disease status

For each file containg predicted mtDNA adundance for each mtDNA-encoded gene in each tissue, test predicted values against quantitative traits and disease status using linear models:

```
for i in *rdat.csv; do stub=$(echo $i | sed 's/.rdat.csv//');
perl test_only_local_full.pl ${stub}.rdat.csv;
done
```
