# prerequisite in path: featureCounts, STAR, mixcr
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

# bam file ("index:bam,bam_file") or fastq files ("index:fastq1,fastq2")
# gtf file
# output folder
# number of threads to use
my ($bam_file,$gtf,$output_folder,$thread)=@ARGV;
my ($path,$index,$fastq);
my $parameters=" --primary -O -t exon -g transcript_id -s 0 -T ".$thread." --largestOverlap --minOverlap 3 --ignoreDup -p -P -B -C";
my $compress=" --readFilesCommand zcat";

#################  clean data  ###############################

$path=abs_path($0);
$path=~s/expression\.pl//;

$bam_file=~/^(.*)\:(.*)$/;
$index=$1;
$fastq=$2;

# if giving bam files
if ($fastq=~s/^bam,//)
{ 
  system("perl ".$path."/bam2fastq.pl ".$fastq." ".$output_folder." ".$thread);
  $fastq=$output_folder."/fastq1.fastq ".$output_folder."/fastq2.fastq";
  $compress="";
}else # if giving fastq file(s)
{
  $fastq=~s/,/ /;
}

#################  expression analyses  ######################

# STAR alignment
system("STAR --runThreadN ".$thread." --genomeDir ".$index." --readFilesIn ".$fastq." --outFileNamePrefix ".$output_folder.
  "/ --outSAMtype BAM Unsorted".$compress);
$bam_file=$output_folder."/Aligned.out.bam";

# featureCounts
unless (-f $gtf) {die "Error: Gtf annotation file doesn't exists!\n";}
unless (-f $bam_file) {die "Error: RNA-Seq bam file doesn't exists!\n";}

system("cd ".$output_folder." && sambamba view -h -t ".$thread." ".$bam_file." > aligned.sam");
system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/transcript.featureCounts aligned.sam ".$parameters);
system("cd ".$output_folder." && featureCounts -a ".$gtf." -o ".$output_folder."/exon.featureCounts aligned.sam -f ".$parameters);
system("rm -f ".$output_folder."/aligned.sam");

unless (-f $output_folder."/transcript.featureCounts") {die "Error: featureCounts failed!\n";}
unless (-f $output_folder."/exon.featureCounts") {die "Error: featureCounts failed!\n";}

# cleanup
unlink($output_folder."/transcript.featureCounts.summary");
unlink($output_folder."/exon.featureCounts.summary");
unlink($output_folder."/Log.out");
unlink($output_folder."/Log.progress.out");
unlink($output_folder."/SJ.out.tab");
unlink($output_folder."/Log.final.out");
unlink($bam_file);

###################  T cell repertoire analyses  #################

system("mixcr align -s hs -f -p rna-seq -t ".$thread." -OallowPartialAlignments=true ".$fastq." ".
  $output_folder."/alignments.vdjca > ".$output_folder."/mixcr1.txt");
system("mixcr assemblePartial -f ".$output_folder."/alignments.vdjca ".$output_folder.
  "/alignment_contigs.vdjca > ".$output_folder."/mixcr21.txt");
system("mixcr assemblePartial -f ".$output_folder."/alignment_contigs.vdjca ".$output_folder.
  "/alignment_contigs2.vdjca > ".$output_folder."/mixcr22.txt");
system("mixcr extend -f ".$output_folder."/alignment_contigs2.vdjca ".$output_folder.
  "/alignmentsRescued_2_extended.vdjca > ".$output_folder."/mixcr3.txt");
system("mixcr assemble -ObadQualityThreshold=15 -OaddReadsCountOnClustering=true -f -t ".
  $thread." ".$output_folder."/alignmentsRescued_2_extended.vdjca ".
  $output_folder."/clones.clns > ".$output_folder."/mixcr4.txt");
system("mixcr exportClones -f ".$output_folder."/clones.clns ".$output_folder.
  "/clones.txt > ".$output_folder."/mixcr5.txt");

unlink($output_folder."/alignments.vdjca");
unlink($output_folder."/alignment_contigs.vdjca");
unlink($output_folder."/alignment_contigs2.vdjca");
unlink($output_folder."/alignmentsRescued_2_extended.vdjca");
unlink($output_folder."/clones.clns");
unlink($output_folder."/mixcr1.txt");
unlink($output_folder."/mixcr21.txt");
unlink($output_folder."/mixcr22.txt");
unlink($output_folder."/mixcr3.txt");
unlink($output_folder."/mixcr4.txt");
unlink($output_folder."/mixcr5.txt");

