###########  nucleus-specific sbatch wrapper for calling somatic mutations  #############
# for other job submission system, you should be able to easily change "sbatch" to the appropriate command
#
# input format:
# jobs: the batch job design file, it has 6 columns separated by \t, they correspond to the $somatic, $expression_somatic, 
#       $output, $fastq1, $fastq2, and $exp_bam input variables of detect_neoantigen.pl
#       (1) Commented lines ("#" at the front) are skipped
#       (2) if providing RNA-Seq fastq file(s) in lieu of $exp_bam, 
#           one STAR aligner index is already provided in the hg38 folder
# example: the demo job submission shell script. A default one is in this folder
# max_normal_cutoff, min_tumor_cutoff, build, gtf, mhc_i, mhc_ii, percentile_cutoff, rpkm_cutoff, thread, max_mutations: follow detect_neoantigen.pl
# n: bundle $n somatic calling jobs into one submission
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($jobs,$example,$max_normal_cutoff,$min_tumor_cutoff,$build,$gtf,$mhc_i,$mhc_ii,$percentile_cutoff, 
  $rpkm_cutoff,$thread,$max_mutations,$n)=@ARGV;
my ($line,$line1,@items,$i,$job);

my $path=abs_path($0);
$path=~s/job_detect_neoantigen\.pl//;

open(JOB,$jobs) or die "Cannot find the design file!\n";
$i=0;

while ($line=<JOB>)
{
  $line=~s/(\r|\n)//g;
  if ($line eq "" || $line=~/^#/) {next;}

  if ($i++ % $n==0)
  {
    $job="neoantigen_".$i.".sh";
    open(SCRIPT,">".$job) or die "Cannot write to the shell submission script!\n";

    # write header
    open(HEADER,$example) or die "Cannot find the example shell script!\n";
    while ($line1=<HEADER>)
    {
      if ($line1=~/JOBSTART/) {last;}
      print SCRIPT $line1;
    }
    close(HEADER);
  }

  # write submission job
  @items=split("\t",$line);
  print SCRIPT "perl ".$path."/detect_neoantigen.pl ".$items[0]." ".$items[1]." ".$max_normal_cutoff." ".
    $min_tumor_cutoff." ".$build." ".abs_path($items[2])." ".$items[3]." ".$items[4]." ".$items[5]." ".$gtf.
    " ".$mhc_i." ".$mhc_ii." ".$percentile_cutoff." ".$rpkm_cutoff." ".$thread." ".$max_mutations."\n";

  if ($i % $n==0)
  {
    close(SCRIPT); 
    system("sbatch ".$job);
    unlink($job);
    sleep(1);
  }
}

close(JOB);

if ($i % $n!=0)
{
  close(SCRIPT);
  system("sbatch ".$job);
  unlink($job);
}

# perl /home2/twang6/software/immune/neoantigen/job_detect_neoantigen.pl \
# design.txt \
# /project/bioinformatics/Xiao_lab/shared/neoantigen/code/neoantigen/example/example.sh \
# 0.02 0.05 hg38 \
# /home2/twang6/data/genomes/hg38/hg38_genes.gtf \
# /project/bioinformatics/Xiao_lab/shared/neoantigen/code/mhc_i \
# /project/bioinformatics/Xiao_lab/shared/neoantigen/code/mhc_ii \
# 3 1 32 50000 2

