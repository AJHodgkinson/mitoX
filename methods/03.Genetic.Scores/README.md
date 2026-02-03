# Generation of Genetic Scores for mtDNA-encoded transcript abundance

Loop through tissues and genes, using mitoimpute (https://github.com/AJHodgkinson/mitoimpute) to create score models.  The following loop was used to generate main European models (printing to shell files to run in parrallel), but parameters can be modified to run different options: 

```
for i in "EUR.AdiposeSubcutaneous.BASIC" "EUR.AdiposeVisceralOmentum.BASIC" "EUR.AdrenalGland.BASIC" "EUR.ArteryAorta.BASIC" "EUR.ArteryCoronary.BASIC" "EUR.ArteryTibial.BASIC" "EUR.BrainAmygdala.BASIC" "EUR.BrainAnteriorcingulatecortexBA24.BASIC" "EUR.BrainCaudatebasalganglia.BASIC" "EUR.BrainCerebellarHemisphere.BASIC" "EUR.BrainCerebellum.BASIC" "EUR.BrainCortex.BASIC" "EUR.BrainFrontalCortexBA9.BASIC" "EUR.BrainHippocampus.BASIC" "EUR.BrainHypothalamus.BASIC" "EUR.BrainNucleusaccumbensbasalganglia.BASIC" "EUR.BrainPutamenbasalganglia.BASIC" "EUR.BrainSpinalcordcervicalc1.BASIC" "EUR.BrainSubstantianigra.BASIC" "EUR.BreastMammaryTissue.BASIC" "EUR.CellsCulturedfibroblasts.BASIC" "EUR.CellsEBVtransformedlymphocytes.BASIC" "EUR.ColonSigmoid.BASIC" "EUR.ColonTransverse.BASIC" "EUR.EsophagusGastroesophagealJunction.BASIC" "EUR.EsophagusMucosa.BASIC" "EUR.EsophagusMuscularis.BASIC" "EUR.HeartAtrialAppendage.BASIC" "EUR.HeartLeftVentricle.BASIC" "EUR.KidneyCortex.BASIC" "EUR.Liver.BASIC" "EUR.Lung.BASIC" "EUR.MinorSalivaryGland.BASIC" "EUR.MuscleSkeletal.BASIC" "EUR.NerveTibial.BASIC" "EUR.Ovary.BASIC" "EUR.Pancreas.BASIC" "EUR.Pituitary.BASIC" "EUR.Prostate.BASIC" "EUR.SkinNotSunExposedSuprapubic.BASIC" "EUR.SkinSunExposedLowerleg.BASIC" "EUR.SmallIntestineTerminalIleum.BASIC" "EUR.Spleen.BASIC" "EUR.Stomach.BASIC" "EUR.Testis.BASIC" "EUR.Thyroid.BASIC" "EUR.Uterus.BASIC" "EUR.Vagina.BASIC" "EUR.WholeBlood.BASIC";
do for j in "ENSG00000198695.2_med_mtnuc" "ENSG00000198712.1_med_mtnuc" "ENSG00000198727.2_med_mtnuc" "ENSG00000198763.3_med_mtnuc" "ENSG00000198786.2_med_mtnuc" "ENSG00000198804.2_med_mtnuc" "ENSG00000198840.2_med_mtnuc" "ENSG00000198886.2_med_mtnuc" "ENSG00000198888.2_med_mtnuc" "ENSG00000198899.2_med_mtnuc" "ENSG00000198938.2_med_mtnuc" "ENSG00000210082.2_med_mtnuc" "ENSG00000211459.2_med_mtnuc" "ENSG00000212907.2_med_mtnuc" "ENSG00000228253.1_med_mtnuc";
do for k in "1e-5" "1e-1" "1";
altered=$(echo $i | sed 's/EUR./EUR_/');
nextflow run main.nf --name "GTEx.Models.${i}.${j}.${k}" --bfile "GTEx_${altered}_maf5miss1hwe" --eqtl "GTEx_${altered}.${j}.assoc.linear.gz" --pvalue "${k}" --pheno "${j}" --phenofile "plink_pheno_GTEx_${altered}.txt" --matchBim "GTEx_${altered}_maf5miss1hwe.bim" --useMito "1" --useCovars "1" --liftmatch "0" --liftcompare "0" --model "lasso" --targetbfile "GTEx_${altered}_maf5miss1hwe" --targetpheno "${j}" --targetphenofile "plink_pheno_GTEx_${altered}.txt" --targetcovars "GPC1,GPC2,GPC3,GPC4,GPC5,PEER1,PEER2,PEER3,PEER4,PEER5,PEER6,PEER7,PEER8,PEER9,PEER10" --useRef "hg20" -profile singularity;
done;
done;
done
```


