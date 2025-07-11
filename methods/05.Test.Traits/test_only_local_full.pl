use strict; 
use Statistics::R;

my $R = Statistics::R->new();
my $make;

my $input=$ARGV[0]; #Name for expression prediction file

my $name=$input; $name=~s/.rdat.csv//;
open (OUT, ">${name}.lm.txt") || die "Unable to open results file to write to: $!\n";

my @disease=(); my @quant=();

if (!($input=~/EUR/)) {
  #Load all data matrices and predicted gene expression:
  $make="disease <- read.delim('/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.self.unrel.txt', header = TRUE, sep = '\t')";
  $R->send($make);
  print "Disease loaded\n";
  $make="quant <- read.delim('/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.quant.unrel.normalised.blood.txt', header = TRUE, sep = '\t')";
  $R->send($make);
  print "Quant loaded\n";
  $make="imputed <- read.csv('/scratch/grp/hodgkinsonlab/shared/GTEx_Data/GTEx_Models/$input', header = TRUE)";
  $R->send($make);
  print "Expression loaded\n";
  $make="covars <- read.delim('/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.covars.unrel.txt', header = TRUE, sep = '\t')";
  $R->send($make);
  print "Covars loaded\n";
  $make="names(imputed)[names(imputed) == \"sample\"] <- \"IID\""; #rename
  $R->send($make);
  
  ###For disease traits first:
  #Get diseases in file:
  open (DISEASE, "/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.self.unrel.txt") || die "Unable to open disease file to read: $!\n";
  while (<DISEASE>) {
    if ($_=~/^FID/) {
      my @array=split;
      for (my $i=2;$i<@array;$i++) {
	my $dd=$array[$i];
	$dd=~s/\-/\./g;
	push @disease, $dd;
      }
      goto endcollect1a;
    }
  }
  
 endcollect1a:
  close (DISEASE);

  ###For QUANT traits first:
  #Get QUANT in file:
  open (QUANT, "/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.quant.unrel.normalised.blood.txt") || die "Unable to open quant file to read: $!\n";
  while (<QUANT>) {
    if ($_=~/^FID/) {
      my @array=split;
      for (my $i=2;$i<@array;$i++) {
	my $dd=$array[$i];
	$dd=~s/\-/\./g;
	push @quant, $dd;
      }
      goto endcollect1b;
    }
  }
  
 endcollect1b:
  close (QUANT);

  print "Features collected\n";
  
}

if ($input=~/EUR/) {
  #Load all data matrices and predicted gene expression:
  $make="disease <- read.delim('/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.self.unrel.eur.txt', header = TRUE, sep = '\t')";
  $R->send($make);
  print "Disease loaded\n";
  $make="quant <- read.delim('/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.quant.unrel.eur.normalised.blood.txt', header = TRUE, sep = '\t')";
  $R->send($make);
  print "Quant loaded\n";
  $make="imputed <- read.csv('/scratch/grp/hodgkinsonlab/shared/GTEx_Data/GTEx_Models/$input', header = TRUE)";
  $R->send($make);
  print "Expression loaded\n";
  $make="covars <- read.delim('/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.covars.unrel.eur.txt', header = TRUE, sep = '\t')";
  $R->send($make);
  print "Covars loaded\n";
  $make="names(imputed)[names(imputed) == \"sample\"] <- \"IID\""; #rename
  $R->send($make);
  
  ###For disease traits first:
  #Get diseases in file:
  open (DISEASE1, "/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.self.unrel.eur.txt") || die "Unable to open disease file to read: $!\n";
  while (<DISEASE1>) {
    if ($_=~/^FID/) {
      my @array=split;
      for (my $i=2;$i<@array;$i++) {
	my $dd=$array[$i];
	$dd=~s/\-/\./g;
	push @disease, $dd;
      }
      goto endcollect2a;
    }
  }
  
 endcollect2a:
  close (DISEASE1);

  ###For QUANT traits first:
  #Get QUANT in file:
  open (QUANT1, "/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/phenotypes/master.file.quant.unrel.eur.normalised.blood.txt") || die "Unable to open quant file to read: $!\n";
  while (<QUANT1>) {
    if ($_=~/^FID/) {
      my @array=split;
      for (my $i=2;$i<@array;$i++) {
	my $dd=$array[$i];
	$dd=~s/\-/\./g;
	push @quant, $dd;
      }
      goto endcollect2b;
    }
  }
  
 endcollect2b:
  close (QUANT1);
  print "Features collected\n";
}

#Add age squared to covars:
$make="covars\$Agesq <- covars\$Age^2";
$R->send($make);



#Make combined data frame:
$make="merged_inter <- merge(disease, imputed, by = \"IID\")";
$R->send($make);
$make="merged_disease <- merge(merged_inter, covars, by = \"IID\")";
$R->send($make);

print "Disease data merged\n";

#Loop through diseases and test:
foreach my $disease (@disease) {
  print "$disease\n";
  $make="results<-lm(merged_disease\$imputed ~ merged_disease\$$disease+merged_disease\$GPC1+merged_disease\$GPC2+merged_disease\$GPC3+merged_disease\$GPC4+merged_disease\$GPC5+merged_disease\$GPC6+merged_disease\$GPC7+merged_disease\$GPC8+merged_disease\$GPC9+merged_disease\$GPC10+merged_disease\$Age+merged_disease\$GenoBatch+merged_disease\$Sex)";
  $R->send($make);
  $make="pval=summary(results)\$coefficients[2,4]";
  $R->send($make);
  my $pval = "NA";
  $pval=$R->get('pval');
  $make="est=summary(results)\$coefficients[2,1]";
  $R->send($make);
  my $est = "NA";
  $est=$R->get('est');
  $make="cases=sum(merged_disease\$$disease==1)";
  $R->send($make);
  my $cases = "NA";
  $cases=$R->get('cases');
  $make="total=nobs(results)";
  $R->send($make);
  my $total = "NA";
  $total=$R->get('total');
  print OUT "$name\t$disease\t$cases\t$total\t$est\t$pval\n";
}

print "Disease data tested\n";

my %set1 = ("Height"=>0,"Weight1"=>0,"Waist_circumference"=>0,"Hip_circumference"=>0,"Waist-to-hip_ratio"=>0,"Forced_expiratory_volume_FEV1"=>0,"Forced_expiratory_vital_capacity_FVC"=>0,"Forced_expiratory_volume_persantage_FEV1PP"=>0,"Airway_function"=>0,"Peak_Expiratory_Flow_PEF"=>0,"Pulse_wave_arterial_stiffness_index"=>0,"Systolic_blood_pressure"=>0,"Diastolic_blood_pressure"=>0,"Pulse_pressure"=>0,"Maximum_heart_rate_during_fitness_test"=>0,"Maximum_workload_during_fitness_test"=>0,"ECG_heart_rate"=>0,"Cardiac_output"=>0,"Intraocular_pressure"=>0,"LogMAR_right_eye"=>0,"LogMAR_left_eye"=>0,"Children"=>0,"Birth_weight"=>0,"Overall_acceleration_average"=>0,"No_wear_time_bias_adjusted_average_acceleration"=>0,"Hand_grip_strength_left"=>0,"Hand_grip_strength_right"=>0); #array, sex, age, age squared, 10 nucPCs
my %set2 = ("Weight2"=>0,"Body_mass_index_BMI"=>0,"Leg_fat_percentage"=>0,"Leg_fat_mass"=>0,"Leg_fat-free_mass"=>0,"Leg_predicted_mass"=>0,"Arm_fat_percentage"=>0,"Arm_fat_mass"=>0,"Arm_fat-free_mass"=>0,"Arm_predicted_mass"=>0,"Trunk_fat_percentage_TFP"=>0,"Trunk_fat_mass_TFM"=>0,"Trunk_fat_free_mass_TFFM"=>0,"Trunk_predicted_mass_TPM"=>0,"Body_fat_percentage_BFP"=>0,"Whole_body_fat_mass_WBFM"=>0,"Whole_body_fat_free_mass_WBFFM"=>0,"Impedance_of_whole_body_IWB"=>0,"Impedance_of_leg"=>0,"Impedance_of_arm"=>0,"Fat_mass_index"=>0,"Fat_free_mass_index"=>0,"Alanine_aminotransferase_ALT"=>0,"Albumin"=>0,"Alkaline_phosphatase_ALP"=>0,"Apolipoprotein_A1"=>0,"Apolipoprotein_B"=>0,"Aspartate_aminotransferase_AST"=>0,"AST_ALT_ratio"=>0,"C_reactive_protein_CRP"=>0,"Calcium"=>0,"Cholesterol_total"=>0,"Creatinine"=>0,"Cystatin_C"=>0,"Direct_bilirubin"=>0,"Gamma_glutamyltransferase_GGT"=>0,"Glucose"=>0,"Glycated_haemoglobin_HbA1c"=>0,"HDL_cholesterol_HDL_C"=>0,"IGF_1"=>0,"LDL_direct"=>0,"Lipoprotein_A"=>0,"Non-HDL"=>0,"Oestradiol"=>0,"Phosphate"=>0,"Rheumatoid_factor"=>0,"Sex_hormone_binding_globulin_SHBG"=>0,"Total_bilirubin"=>0,"Testosterone"=>0,"Total_protein"=>0,"Triglycerides"=>0,"Urate"=>0,"Urea"=>0,"Vitamin_D"=>0,"Estimated_Glomerular_Filtration_Rate_creatinine"=>0,"Estimated_Glomerular_Filtration_Rate_Cystatin"=>0,"Estimated_Glomerular_Filtration_Rate_creatinine_cystatin"=>0,"Microalbumin_in_urine"=>0,"Creatinine_enzymatic_in_urine"=>0,"Potassium_in_urine"=>0,"Sodium_in_urine"=>0,"Basal_Metabolic_Rate"=>0); #array, sex, 10 nucPCs
my %set3=("Platelet_count"=>0,"White_blood_cell_count"=>0,"Red_blood_cell_count"=>0,"Hemoglobin_concentration"=>0,"Hematocrit"=>0,"Mean_corpuscular_volume"=>0,"Mean_corpuscular_hemoglobin"=>0,"Mean_corpuscular_hemoglobin_concentration"=>0,"Red_cell_distribution_width"=>0,"Plateletcrit"=>0,"Mean_platelet_volume"=>0,"Platelet_distribution_width"=>0,"Lymphocyte_count"=>0,"Monocyte_count"=>0,"Neutrophil_count"=>0,"Eosinophil_count"=>0,"Basophil_count"=>0,"Lymphocyte_percentage_of_white_cells"=>0,"Monocyte_percentage_of_white_cells"=>0,"Neutrophil_percentage_of_white_cells"=>0,"Eosinophil_percentage_of_white_cells"=>0,"Basophil_percentage_of_white_cells"=>0,"Reticulocyte_fraction_of_red_cells"=>0,"Reticulocyte_count"=>0,"Immature_fraction_of_reticulocytes"=>0,"High_light_scatter_reticulocyte_percentage_of_red_cells"=>0,"High_light_scatter_reticulocyte_count"=>0,"Granulocyte_count"=>0,"Neutrophil_percentage_of_granulocytes"=>0,"Basophil_percentage_of_granulocytes"=>0,"Eosinophil_percentage_of_granulocytes"=>0,"Myeloid_white_cell_count"=>0,"Granulocyte_percentage_of_myeloid_white_cells"=>0,"Age_mother"=>0,"Age_father"=>0,"Mean_age"=>0,"AgeZ"=>0,"Reproductive_span"=>0,"Age_at_menopause"=>0,"Age_at_menarche"=>0); #array, 10 nucPCs

  
#Make combined data frame:
$make="merged_inter1 <- merge(quant, imputed, by = \"IID\")";
$R->send($make);
$make="merged_quant <- merge(merged_inter1, covars, by = \"IID\")";
$R->send($make);

print "Quant data merged\n";

#Loop through quant and test:
foreach my $trait (@quant) {
  print "$trait\n";
  $make="results<-lm(merged_quant\$imputed ~ as.numeric(merged_quant\$$trait)+merged_quant\$GPC1+merged_quant\$GPC2+merged_quant\$GPC3+merged_quant\$GPC4+merged_quant\$GPC5+merged_quant\$GPC6+merged_quant\$GPC7+merged_quant\$GPC8+merged_quant\$GPC9+merged_quant\$GPC10+merged_quant\$Age+merged_quant\$Agesq+merged_quant\$GenoBatch+merged_quant\$Sex)" if (exists $set1{$trait});
  $make="results<-lm(merged_quant\$imputed ~ as.numeric(merged_quant\$$trait)+merged_quant\$GPC1+merged_quant\$GPC2+merged_quant\$GPC3+merged_quant\$GPC4+merged_quant\$GPC5+merged_quant\$GPC6+merged_quant\$GPC7+merged_quant\$GPC8+merged_quant\$GPC9+merged_quant\$GPC10+merged_quant\$GenoBatch+merged_quant\$Sex)" if (exists $set2{$trait});
  $make="results<-lm(merged_quant\$imputed ~ as.numeric(merged_quant\$$trait)+merged_quant\$GPC1+merged_quant\$GPC2+merged_quant\$GPC3+merged_quant\$GPC4+merged_quant\$GPC5+merged_quant\$GPC6+merged_quant\$GPC7+merged_quant\$GPC8+merged_quant\$GPC9+merged_quant\$GPC10+merged_quant\$GenoBatch)" if (exists $set3{$trait}); 
  $R->send($make);
  $make="pval=summary(results)\$coefficients[2,4]";
  $R->send($make);
  my $pval = "NA";
  $pval=$R->get('pval');
  $make="est=summary(results)\$coefficients[2,1]";
  $R->send($make);
  my $est = "NA";
  $est=$R->get('est');
  $make="total=nobs(results)";
  $R->send($make);
  my $total = "NA";
  $total=$R->get('total');
  print OUT "$name\t$trait\t$total\t$total\t$est\t$pval\n";
}

close (OUT);

print "Quant data tested\n";
