use strict;

#Generate a file-manifest for samples in the study, and then extract each file name and object_id for download with gen3-client
open (FILE, "file-manifest.json") || die "Unable to open file to read: $!\n";

my $rna; my $dg; my $md5;

while (<FILE>) {
  my @array=split;
  if ($_=~/md5sum/) { #Get md5sum to check download against. The file name is in the second column for lines containing "md5sum".
    $md5=$array[1];
    $md5=~s/\"//g;
    $md5=~s/\,//g;
  }
  if ($_=~/file_name/) { #Collect RNA sequencing filename. The file name is in the second column for lines containing "file_name".
    $rna=$array[1];
    $rna=~s/\"//g;
    $rna=~s/\,//g;
  }
  if ($_=~/object_id/) { #Collect object ID. The object name is in the second column for lines containing "object_id".
    $dg=$array[1];
    $dg=~s/\"//g;
    $dg=~s/\,//g;
    
    my @out=split(/\./,$rna);
    my $stub=$out[0]; #Collect the root name of the sample bu splitting on full stop and taking the first element
    my $outfile=$stub.".sh"; #Create out shell file for each sample
    my $check=$stub.".md5"; #Create the md5sum filename
    my $complete=$stub.".MT.SEfile.txt"; #Create the output file name
    
    if (!( -e $complete)) { #If the proposed output file does not already exist:

      #Open the md5 output file
      open (CHE, ">$check") || die "Unable to open file to read: $\n";
      print CHE "$md5\t$rna\n"; #write the md5sum and rna filename to md5 output file
      close (CHE);
      
      open (OUT, ">$outfile") || die "Unable to open file to read: $\n"; #open the shell file to write to:

      #Add commands to the shell file to use bash, to output the appropriate report file and to change to tthe he correct directory
      print OUT "#!/bin/bash -l\n#SBATCH --output=/scratch/grp/hodgkinsonlab/shared/GTEx_Data/RNAseq/$stub.report\n";
      print OUT "source ~/.bashrc\n";
      print OUT "cd /scratch/grp/hodgkinsonlab/shared/GTEx_Data/RNAseq\n";

      #Instructions to download each file, using gen3 client and the required file details
      print OUT "yes | gen3-client download-single --profile=ahodgkinson --guid=$dg --download-path=/scratch/grp/hodgkinsonlab/shared/GTEx_Data/RNAseq --protocol=s3\n";
      print OUT "md5sum -c $check\n";
      print OUT "rm $check\n";

      #Load samtools for use in parsing the data
      print OUT "module load samtools/1.14-gcc-10.3.0-python-2.7.18\n";

      #Index downloaded BAM file, convert to fastq using samtools.  Original data is BAM file containing all reads (including unmapped reads).
      print OUT "samtools index $rna\n";
      print OUT "samtools sort -@ 10 -n -o $stub.sorted.bam $rna\n"; #sort the file by read
      print OUT "samtools fastq $stub.sorted.bam -1 ${stub}_1.fq.gz -2 ${stub}_2.fq.gz\n"; #convert to fastq format
      print OUT "rm $stub.sorted.bam\n";
      print OUT "rm $rna\n";
      print OUT "rm ${rna}.bai\n";

      #print processign commands to shell file: Align with STAR, generate gene count data, filter properly paired, uniquely mapped reads, and create MT only file
      print OUT "module load star/2.7.6a-gcc-9.4.0\n"; #load STAR v2.7.6
      #print star command: use predefined star reference, quantmode to get gene counts, 2-pass mode and defining overlap, permitted mismatch rate and flags to report (see star manual for more details on these options).
      print OUT "STAR --runThreadN 10 --genomeDir /scratch/grp/hodgkinsonlab/shared/GTEx_Data/References/STARGenome --readFilesIn ${stub}_1.fq.gz ${stub}_2.fq.gz --outFileNamePrefix ${stub}. --quantMode GeneCounts --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000  --alignMatesGapMax 1000000 --outSAMattributes NH nM NM MD HI\n";
      print OUT "samtools index ${stub}.Aligned.sortedByCoord.out.bam\n";
      #Keep only properly paired reads with 0x0002 sam flag:
      print OUT "samtools view -@ 10 -f 0x0002 -b -o ${stub}.Aligned.sortedByCoord.out.PP.bam ${stub}.Aligned.sortedByCoord.out.bam\n";
      print OUT "samtools index ${stub}.Aligned.sortedByCoord.out.PP.bam\n";
      #Keep only uniquely mapped reads by filtering on NH:i:1 (one best match):
      print OUT "samtools view -h ${stub}.Aligned.sortedByCoord.out.PP.bam | grep -P \"NH:i:1\\t|^@\" | samtools view -bS - > ${stub}.Aligned.sortedByCoord.out.PP.UM.bam\n";
      print OUT "samtools index ${stub}.Aligned.sortedByCoord.out.PP.UM.bam\n";
      #Keep only data coming from mtDNA
      print OUT "samtools view -bh ${stub}.Aligned.sortedByCoord.out.PP.UM.bam chrM >> ${stub}.Aligned.sortedByCoord.out.PP.UM.MT.bam\n";
      print OUT "samtools index ${stub}.Aligned.sortedByCoord.out.PP.UM.MT.bam\n";

      #Create modification file, cleavage file and QC data
      #Create read start/end file - outputs number of reads starting and ending at each position in the mtDNA genome:
      print OUT "perl /scratch/grp/hodgkinsonlab/shared/GTEx_Data/References/ReadSEExtractor_mito.pl --Bam ${stub}.Aligned.sortedByCoord.out.PP.UM.MT.bam --RefFasta /scratch/grp/hodgkinsonlab/shared/GTEx_Data/References/GRCh38.primary_assembly.genome.fa --Out ${stub}\n";
      #Create nucleotide count file for mtDNA - outputs number of each nucleotide (forward and reverse) at each position in the mtDNA genome:
      print OUT "perl /scratch/grp/hodgkinsonlab/shared/GTEx_Data/References/pileupAlleleExtractor_mito.pl --Bam ${stub}.Aligned.sortedByCoord.out.PP.UM.MT.bam --MinQ 23 --RefFasta /scratch/grp/hodgkinsonlab/shared/GTEx_Data/References/GRCh38.primary_assembly.genome.fa --Out ${stub}\n";
      #Run RNAseQC on BAM file to get quality metrics for filtering:
      print OUT "/scratch/grp/hodgkinsonlab/Programs/rnaseqc.v2.4.2.linux /scratch/grp/hodgkinsonlab/shared/GTEx_Data/References/gencode.v42.primary_assembly.genes.gtf ${stub}.Aligned.sortedByCoord.out.PP.bam RNAseQC\n";
      #Clean up data by removing all intermediate files:
      print OUT "rm -r ${stub}._STARpass1\n";
      print OUT "rm -r ${stub}._STARgenome\n";
      print OUT "rm ${stub}.Log.progress.out\n";
      print OUT "rm ${stub}.Log.out\n";
      print OUT "rm ${stub}.Log.final.out\n";
      print OUT "#rm ${stub}.Aligned.sortedByCoord.out.bam\n";
      print OUT "#rm ${stub}.Aligned.sortedByCoord.out.bam.bai\n";
      print OUT "rm ${stub}.Aligned.sortedByCoord.out.PP.bam\n";
      print OUT "rm ${stub}.Aligned.sortedByCoord.out.PP.bam.bai\n";
      print OUT "rm ${stub}.Aligned.sortedByCoord.out.PP.UM.bam\n";
      print OUT "rm ${stub}.Aligned.sortedByCoord.out.PP.UM.bam.bai\n";
      print OUT "rm ${stub}.SJ.out.tab\n";
      print OUT "rm ${stub}_1.fq.gz\n";
      print OUT "rm ${stub}_2.fq.gz\n";
    }
  }
}
  
close (FILE);
  
