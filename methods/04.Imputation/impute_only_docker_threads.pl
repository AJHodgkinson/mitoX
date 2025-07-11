use strict; 
use Statistics::R;
use threads;
use Thread::Queue;

my $R = Statistics::R->new();
my $make;

my $covars="NA";

my $name=$ARGV[0]; #Name for output files
my $rweights=$ARGV[1]; #R results file from Fusion
my $model=$ARGV[2]; #Select the model for testing
my $plink_stub=$ARGV[3]; #Plink binary stub for target population
my $lift=$ARGV[4]; #0 if not required, hg19tohg20 or hg20tohg19 (fusion data genome on the left, target data genome on the right)
my $selected=$ARGV[5]; #Bim file from selected SNPs so number can be recorded
my $num_threads = $ARGV[6]; #Number of threads to use

my $queue = Thread::Queue->new();

my $report_header="Name";
my $report="$name";

#Open BIM and count selected SNPs:
my $starting=0;
open (STARTING, "$selected") || die "Unable to open weights file to read: $!\n";
while (<STARTING>) {
  $starting++;
}
close (STARTING);

$report_header.="\tStartingVariants";
$report.="\t$starting";

#Open Fusion model output
$make="data<-load(\"$rweights\")";
$R->send($make);

##Store models tested an accompanying statistics
my @models=(); my @rsq=(); my @pval=(); 
my $performance = $R->get('cv.performance');
my $state=0;
foreach my $element (@{$performance}) {
  push @models, $element if (($state==0)&&(!($element=~/rsq|pval/)));
  push @rsq, $element if (($state==1)&&(!($element=~/rsq|pval/)));
  push @pval, $element if (($state==2)&&(!($element=~/rsq|pval/)));
  $state++ if ($element=~/rsq|pval/);
}

##Calculate position of selected model for downstream marker selection
my $position;
for (my $i=0;$i<@models;$i++) {
  $position=$i if ($models[$i] eq $model);
}

my $start=@models;

#Store weights and SNPs for selected model
$make="write.table(wgt.matrix, file = \"${name}.table.weights.txt\", sep = \"\t\")";
$R->send($make);
my @weights=();
open (WEIGHTS, "${name}.table.weights.txt") || die "Unable to open weights file to read: $!\n";
while (<WEIGHTS>) {
  chomp;
  my @array=split;
  foreach my $item (@array) {
    $item=~s/\"//g;
    push @weights, $item;
  }
}
close (WEIGHTS);

my $snp_count=0;

my %weights=();
for (my $i=$start;$i<@weights;$i=$i+($start+1)) {
  $weights{$weights[$i]}=$weights[$i+1+$position];
  $snp_count++ if ($weights[$i+1+$position] != 0);
}


#Store SNP data for selected model
$make="write.table(snps, file = \"${name}.table.snps.txt\", sep = \"\t\", row.names = FALSE)";
$R->send($make);
my @snps=();
open (SNPS, "${name}.table.snps.txt") || die "Unable to open SNPs file to read: $!\n";
while (<SNPS>) {
  chomp;
  my @array=split;
  foreach my $item (@array) {
    $item=~s/\"//g;
    push @snps, $item;
  }
}
close (SNPS);

my %coordinates=(); my %weight_allele=(); my %snp_alleles=(); my %collect=();
#print "SNPS:\n";
for (my $i=6;$i<@snps;$i=($i+6)) {
  my $chr=$snps[$i]; $chr=~s/MT/26/; $chr=~s/M/26/;
  $coordinates{$snps[$i+1]}=$chr."_".$snps[$i+3];
  $weight_allele{$snps[$i+1]}=$snps[$i+4];
  $snp_alleles{$snps[$i+1]}=$snps[$i+4]."_".$snps[$i+5];
  my $tag=$chr."_".$snps[$i+3]."_".$snps[$i+4]."_".$snps[$i+5];
  $collect{$tag}=$snps[$i+1];
}
$report_header.="\tVariants";
$report.="\t$snp_count";

#Report CV statistics
my $hsq = $R->get('hsq');
my $middle="NA";
$middle=1 if ($hsq eq "1");

if ($hsq ne "1") {
  my @hsq=@{$hsq};
  my $upper=$hsq[0]; my $lower=$hsq[1]; $middle=$lower+(($upper-$lower)/2);
}

$report_header.="\tHerit";
$report.="\t$middle";

my $hsqpv = $R->get('hsq.pv');
my $tot = $R->get('N.tot');

$report_header.="\tHerit_P\tIndividuals";
$report.="\t$hsqpv\t$tot";

for (my $i=0;$i<@models;$i++) {
  $report_header.="\tModel\t5F-Rsq\t5F-P" if ($models[$i] eq $model);
  $report.="\t$models[$i]\t$rsq[$i]\t$pval[$i]" if ($models[$i] eq $model);
}

#If liftover is required, create BED file for selected SNPs, run liftover, and store SNPs for later comparison
my %liftover=();

open (LIFT, ">$name.lift.bed") || die "Unable to open lift file to write to: $!\n";
for (my $i=6;$i<@snps;$i=($i+6)) {
  my $chr=$snps[$i]; $chr=~s/MT/26/; $chr=~s/M/26/;
  my $pos=$snps[$i+3]; my $pos1=$pos+1;
  print LIFT "chr$chr\t$pos\t$pos1\t${chr}_$pos\n"; #Add 'chr' so liftover will work
}

close (LIFT);

system `liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg38ToHg19.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg19tohg20/);
system `liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg19ToHg38.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg20tohg19/);
system `grep chr26 $name.lift.bed >> $name.lifted.bed` if ($lift=~/hg/); #Add Mito chromosome, as this will not have been lifted over (but co-ordinates are the same between hg19 and hg20)

system `cp $name.lift.bed $name.lifted.bed` if (!($lift=~/hg/)); #If no liftover reuqured, just copy the file
system `touch $name.not_lifted.bed` if (!($lift=~/hg/)); #If no liftover reuqured, just copy the file

open (LIFTED, "$name.lifted.bed") || die "Unable to open $name.lifted.bed to read: $!\n";
while (<LIFTED>) {
  my @array=split;
  my $tag=$array[0]."_".$array[1]; $tag=~s/chr//;
  $liftover{$tag}=$array[3];
}
close (LIFTED);


#####Copy over target SNP file and convert SNP names to be compatible with Fusion output model

open (BIMFULL, "/scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/imputed/ukb.imp.maf1.miss5.hwe.MT.removed.bim") || die "Unable to open ${plink_stub}.bim to read: $!\n";
open (BIMHOLD, ">$name.full.bim") || die "Unable to write to $name.full.bim to read: $!\n";

#my @markers=();
while (<BIMFULL>) {
  chomp;
  my @array=split;
  my $chr=$array[0]; $chr=~s/chr//; $chr=~s/MT/26/; $chr=~s/M/26/;
  print BIMHOLD "$chr\t${chr}_$array[3]_$array[5]_$array[4]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\n";
  my $marker="${chr}_$array[3]_$array[5]_$array[4]";
  #push @markers, $marker;
}

close (BIMFULL);
close (BIMHOLD);

#Loop through BIM file for taregt data and collect SNP information for SNPs in model, and also collect weight and allele information from these sites:
open (BIM, "$name.full.bim") || die "Unable to open target BIM file to read: $!\n";
open (COLLECT, ">$name.snps.txt") || die "Unable to write to SNP file to read: $!\n";
my @feature_weights=(); my @weight_allele=();

while (<BIM>) {
  my @array=split;
  my $chrpos=$array[0]."_".$array[3];
  if (!($lift=~/hg/)) {
    my $tag1=$chrpos."_".$array[4]."_".$array[5];
    my $tag2=$chrpos."_".$array[5]."_".$array[4];
    if (exists $collect{$tag1}) {
      if ($weights{$collect{$tag1}} != 0) {
	print COLLECT "$array[1]\n";
	push @feature_weights, $weights{$collect{$tag1}} if (exists $weights{$collect{$tag1}});
	push @weight_allele, $weight_allele{$collect{$tag1}} if (exists $weight_allele{$collect{$tag1}});
      }
    }
    if (exists $collect{$tag2}) {
      if ($weights{$collect{$tag2}} != 0) {
	print COLLECT "$array[1]\n";
	push  @feature_weights, $weights{$collect{$tag2}} if (exists $weights{$collect{$tag2}});
	push @weight_allele, $weight_allele{$collect{$tag2}} if (exists $weight_allele{$collect{$tag2}});
      }
    }
  }
  if ($lift=~/hg/) { #Use liftover conversion if required
    if (exists $liftover{$chrpos}) {
      my $chrpos1=$liftover{$chrpos};
      my $tag1=$chrpos1."_".$array[4]."_".$array[5];
      my $tag2=$chrpos1."_".$array[5]."_".$array[4];
      if (exists $collect{$tag1}) {
	if ($weights{$collect{$tag1}} != 0) {
	  print COLLECT "$array[1]\n";
	  push  @feature_weights, $weights{$collect{$tag1}} if (exists $weights{$collect{$tag1}});
	  push @weight_allele, $weight_allele{$collect{$tag1}} if (exists $weight_allele{$collect{$tag1}});
	}
      }
      if (exists $collect{$tag2}) {
	if ($weights{$collect{$tag2}} != 0) {
	  print COLLECT "$array[1]\n";
	  push  @feature_weights, $weights{$collect{$tag2}} if (exists $weights{$collect{$tag2}});
	  push @weight_allele, $weight_allele{$collect{$tag2}} if (exists $weight_allele{$collect{$tag2}});
	}
      }
    }
  }  
}
close (BIM);
close (COLLECT);

my $snps_used=@feature_weights;

$report_header.="\tSNPs_Used";
$report.="\t$snps_used";

#####Filter just required SNPs from target file

system `plink2 --bed /scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/imputed/ukb.imp.maf1.miss5.hwe.MT.removed.bed --bim $name.full.bim --fam /scratch/grp/hodgkinsonlab/ahodgkinson/UKBB/UKBB2024/imputed/ukb.imp.maf1.miss5.hwe.MT.removed.fam --extract $name.snps.txt --recode ped --out $name.run`;

#####Impute mitochondria features to target file

eval {

  # Create worker threads
  my @threads;
  for (1 .. $num_threads) {
    push @threads, threads->create(\&process_lines, $name, \@feature_weights, \@weight_allele);
  }
  
  # Open input file for reading
  open my $ped_fh, '<', "$name.run.ped" or die "Unable to open run BIM file to read: $!";
  
  # Enqueue lines for processing
  $queue->enqueue($_) while <$ped_fh>;
  
  # Signal that no more lines will be enqueued
  $queue->end();
  
  # Process the lines and collect results
  my @results = map { $_->join() } @threads;
  
  # Write results to output file
  open my $rdat_fh, '>', "$name.rdat.csv" or die "Unable to write to R file: $!";
  print $rdat_fh "sample,imputed\n";
  print $rdat_fh $_ for @results;
  close $rdat_fh;
  
  # Close the input file handle
  close $ped_fh;
  
  print "Impute collected\n";

  #####Report impute stats
  
  open (RESULTS, ">$name.model.txt") || die "Unable to open results file to write to: $!\n";
  print RESULTS "$report_header\n$report\n";
  close (RESULTS);

};


# Define the processing subroutine
sub process_lines {
    my ($name, $feature_weights_ref, $weight_allele_ref) = @_;
    my $output = '';
    while (defined(my $line = $queue->dequeue())) {
        my @array = split ' ', $line;
        my $impute_sample = 0;
        my $marker = 0;
        for (my $i = 6; $i < @array; $i += 2) {
            $impute_sample += $feature_weights_ref->[$marker] if ($array[$i] eq $weight_allele_ref->[$marker]);
            $impute_sample += $feature_weights_ref->[$marker] if ($array[$i + 1] eq $weight_allele_ref->[$marker]);
            $marker++;
        }
        $output .= "$array[0],$impute_sample\n";
    }
    return $output;
}
