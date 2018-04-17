# prerequisite in path: annovar, humandb of annovar must be in its default position
# somatic: somatic mutation calling file
# normal_freq_cutoff: max VAF in normal sample
# tumor_freq_cutoff: min VAF in tumor sample
# build: hg19 or hg38
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use List::Util qw(max min);

my ($somatic,$normal_freq_cutoff,$tumor_freq_cutoff,$build,$rpkm_cutoff,$max_mutations)=@ARGV;
my ($len,$path,$line,@items,@items1,$header,%aa,$type);
my ($annovar_db,$transcript,$mutation,$class,$extracted);
my %neoantigen_len=("I"=>[8..11],"II"=>[15]); # length of neoantigen for HLA A/B/DRB/DQB

######  filter somatic mutations  ##########

$path=abs_path($0);
$path=~s/predict_neoantigen\.pl//;
system("Rscript ".$path."/filter_somatic.R ".$somatic." ".$tumor_freq_cutoff." ".$normal_freq_cutoff." ".
  $somatic."_filtered ".$rpkm_cutoff." ".$path." ".$max_mutations);
unlink($somatic."_error_no_or_too_many_mutations_after_filtering.txt");
if (! -e $somatic."_filtered") 
{
  system("echo \"\" > ".$somatic."_error_no_or_too_many_mutations_after_filtering.txt");
  exit 1;
}

######  call annotate_variation  ##########

open(AV,"which annotate_variation.pl |") or die "Error: No annovar found!\n";
$annovar_db=<AV>;
close(AV);
$annovar_db=~s/annotate_variation.pl\n/humandb/;
system("annotate_variation.pl -geneanno -dbtype refGene -buildver ".$build." ".$somatic."_filtered ".$annovar_db);

######  expand isoforms of the same gene into multipe lines  ##########

open(FILE_IN,$somatic."_filtered.exonic_variant_function") or die "Error: Didn't find EVF file!\n";
open(FILE_OUT,">".$somatic.".exonic_variant_function_multiple") or die "Error: Cannot write to expanded EVF file!\n";

while ($line=<FILE_IN>)
{
  $line=~s/\n//;
  @items=split("\t",$line);
  $items[1]=~s/ /_/g;
  $items[0].="-".join("-",@items[(1,8..11)]);
  @items1=split(",",$items[2]);
  foreach (0..$#items1)
    {print FILE_OUT join("\t",($items[0]."-".$items1[$_],$items[1],$items1[$_],@items[(3..7)]))."\n";}
}

close(FILE_OUT);
close(FILE_IN);

########  call coding_change  ###########

system("coding_change.pl --includesnp ".$somatic.".exonic_variant_function_multiple ".$annovar_db."/".$build."_refGene.txt ".
  $annovar_db."/".$build."_refGeneMrna.fa > ".$somatic.".coding.txt 2>/dev/null");
system("echo '>WILDTYPE' >> ".$somatic.".coding.txt");
unlink($somatic.".exonic_variant_function_multiple");
unlink($somatic."_filtered.exonic_variant_function");
unlink($somatic."_filtered.log");
unlink($somatic."_filtered.variant_function");

########  extract neoantigens  #########

open(FILE_IN,$somatic.".coding.txt") or die "Error: Didn't find translated proteins!\n";
foreach $class (("I","II"))
{
  foreach $len (@{$neoantigen_len{$class}})
    {unlink($somatic."_".$class."_".$len.".fa");}
}

%aa=("WT"=>"","MU"=>"");
$header="";

while ($line=<FILE_IN>)
{
  if ($line=~/^>/) # this is a header
  {
    if ($line=~/WILDTYPE/) 
    {
      if ($aa{"WT"} ne "" && ($header=~/protein-altering|immediate-stoploss/)) # process the previous AA sequence
      {
        # extract mutation information
        $header=~/(>.*?) /;
        $mutation=$1; 
        $mutation=~s/^>line/>mutation_/;

        foreach $class (("I","II"))
        {
          foreach $len (@{$neoantigen_len{$class}})
          {
            open(FILE_OUT,">>".$somatic."_".$class."_".$len.".fa") or die "Error: Cannot write to fasta files for predicting affinities!\n";
            $extracted=extract_seq($len,$aa{"WT"},$aa{"MU"});
            if (length($extracted)>=$len) {print FILE_OUT $mutation."\n".$extracted."\n";}
            close(FILE_OUT);
          }
        }
      }

      %aa=("WT"=>"","MU"=>"");
      $type="WT";
    }else
    {
      $type="MU";
      $header=$line;
    } 
  }else # append lines to the AA sequence
  {
    $line=~s/\n//;
    $aa{$type}.=$line;
  }
}

close(FILE_IN);

sub extract_seq
{
  my ($overhang,$wt_seq,$mut_seq)=@_;  
  my ($start,$end_wt,$end_mut);

  $wt_seq=~s/\*//;
  $mut_seq=~s/\*//;

  # find starting position
  $start=0;
  while (substr($wt_seq,$start,1) eq substr($mut_seq,$start,1)) 
  {
    $start++;
    if ($start>min(length($wt_seq),length($mut_seq))-1) {last;}
  }

  # find ending position
  $end_wt=length($wt_seq)-1;
  $end_mut=length($mut_seq)-1;
  while (substr($wt_seq,$end_wt,1) eq substr($mut_seq,$end_mut,1) && min($end_wt,$end_mut)>=$start) 
    {$end_wt--;$end_mut--;}
  
  # extract sequence
  $start=max($start-$overhang+1,0);
  $end_mut=min($end_mut+$overhang-1,length($mut_seq)-1);  
  substr($mut_seq,$start,$end_mut-$start+1);
}

