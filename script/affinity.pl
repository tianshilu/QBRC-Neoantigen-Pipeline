# prerequisite: http://www.iedb.org/
# mhc1, mhc2: folders to the iedb mhc1 and mhc2 binding prediction algorithms
# somatic: somatic mutation calling file
# typing: HLA typing results
# max_rank: maximum rank to keep a peptide, recommemded 2
#/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use Scalar::Util qw(looks_like_number);
use Parallel::ForkManager;

my ($mhc1,$mhc2,$somatic,$typing,$max_rank)=@ARGV;
my (@items,$line,$command,$type,$type1,$type2,$class,$file_fa,$file_neoantigen);
my ($pid,$pm,$protein,$path,%protein,%mutations,$i,$len,$len1);
my @dqa1=("","");
my %neoantigen_len=("I"=>[8..11],"II"=>[15]); # length of neoantigen for HLA A/B/DRB/DQB

##########  read existing protein seq  #############

$path=abs_path($0);
$path=~s/affinity\.pl//;

open(PT,$path."/../data/human_protein_seq.fasta") or die "Error: Cannot open the human protein sequence file!\n";

while ($line=<PT>)
{
  if ($line=~/>(.*)\n/)
  {
    $protein=$1;
    $protein{$protein}="";
  }else
  {
    $line=~s/\n//;
    $protein{$protein}.=$line;
  }
}

close(PT);

##########  run iedb analysis  ##############

# extract DQA first
open(TP,$typing) or die "Error: Cannot open the HLA typing file!\n";
while ($line=<TP>)
{
  @items=split("\t",$line);
  if ($items[0] ne "DQA1") {next;}
  @dqa1=($items[1],$items[2]);  
}
close(TP);

# run iedb analysis
open(TP,$typing);
$pm=Parallel::ForkManager->new(8);

while ($line=<TP>)
{
  if ($line=~/DQA1/) {next;} # skip DQA1, only process DQB1
  $pid = $pm -> start and next;

  # parse typing file
  print "Working on ".$line;
  @items=split("\t",$line);
  $type=$items[0];
  $type1=guess_dqa($items[1],$items[2],\@dqa1,0);
  $type2=guess_dqa($items[2],$items[1],\@dqa1,1);

  if ($type eq "A" || $type eq "B" || $type eq "C") 
  {
    $class="I";
    $command=$mhc1."/src/predict_binding.py IEDB_recommended";
  }else 
  {
    $class="II";
    $command=$mhc2."/mhc_II_binding.py NetMHCIIpan";
  }
  
  # run IEDB
  foreach $len (@{$neoantigen_len{$class}})
  {
    # special treatment for class II HLA
    $len1=$len;
    if ($class eq "II") {$len1="";}
    $file_fa=$somatic."_".$class."_".$len.".fa";
    $file_neoantigen=$somatic."_".$type."_".$len."_neoantigen";

    # get mutation information
    $i=1;
    %mutations=("seq_num"=>"mutation\ttype\ttumor_ref\ttumor_var\tnormal_ref\tnormal_var\ttranscript");
    open(FILE_FA,$file_fa) or die "Error: Cannot open the intermediate fasta file prepared for identifying putative neoantigens!\n";
    while ($line=<FILE_FA>) 
    {
      if ($line!~/^>(.*)\n/) {next;}
      $mutations{$i}=$1;
      foreach (1..6) {$mutations{$i}=~s/-/\t/;}
      $i++;
    }
    close(FILE_FA);

    # IEDB
    run_IEDB($command,$type1,$len1,$file_fa,$file_neoantigen."_1.txt",\%mutations,$typing."_".$type);
    run_IEDB($command,$type2,$len1,$file_fa,$file_neoantigen."_2.txt",\%mutations,$typing."_".$type);
  }

  $pm->finish;
}

close(TP);
$pm->wait_all_children;

# clean up
foreach $class (keys %neoantigen_len)
{
  foreach $len (@{$neoantigen_len{$class}})
    {unlink($somatic."_".$class."_".$len.".fa");}
}

system("cat ".$typing."_*invalid | grep \"HLA-\" | sort | uniq > ".$typing.".invalid");
system("rm -f ".$typing."_*invalid");

########  function  ######################

sub run_IEDB
{
  my ($command,$type,$len1,$file_fa,$file_neoantigen,$mutations_ref,$typing_type)=@_;
  my ($allele,$peptide,$rank,@items,$line);
  $allele=0;
  $peptide=5;
  $rank=7;

  open(IEDB,$command." ".$type." ".$len1." ".$file_fa." 2>> ".$typing_type.".invalid |") or die "Error: Cannot execute ".$command."\n";
  open(NAT,">".$file_neoantigen) or die "Error: Cannot write to file: ".$file_neoantigen."!\n";
  while ($line=<IEDB>) # attach mutation information to neoantigen binding strength
  {
    @items=split("\t",$line);
    if ($items[$rank]=~/percentile_rank/) {$items[$rank]="percentile_rank";} # adjust column names
    if ($items[$peptide]=~/core_peptide/) {$items[$peptide]="peptide";}
    $items[$rank]=~s/\n//;
    
    unless ($items[$rank] eq "percentile_rank" || (looks_like_number($items[$rank]) && $items[$rank]<$max_rank)) {next} # not a good binding peptide
    print NAT $mutations_ref->{$items[1]}."\t".match_peptide($items[$peptide],\%protein)."\t".
      join("\t",@items[($allele,$peptide,$rank)])."\n";
  }
  close(IEDB);
  close(NAT);
}

sub match_peptide # see if the putative peptide matches the sequence of any wild type protein
{
  my ($peptide_fa,$protein_ref)=@_;
  my $peptide_fa_regex=qr/$peptide_fa/;
  if ($peptide_fa eq "peptide") {return "wildtype_match";} # header  

  foreach (keys %$protein_ref)
  {
    if ($protein_ref->{$_}=~/$peptide_fa_regex/) {return $_;}
  } 
  return "none_match";
}

sub guess_dqa # match DQA to DQB
{
  my ($type,$type_other,$ref_dqa1,$index)=@_;
  if ($type!~/DQB/) {return $type;} # not DQ allels

  if ($$ref_dqa1[0] ne "") # DQA alleles are typed
  {
    if ($$ref_dqa1[0] eq $$ref_dqa1[1]) # DQA1 alleles are the same
    {
      return $$ref_dqa1[0]."/".$type;
    }else # DQA1 alleles are not the same
    {
      if ($type eq $type_other) # DQB1 alleles are the same
      { 
        return $$ref_dqa1[$index]."/".$type;
      }else # DQB1 alllels are not the same, can only guess
      {
        # http://www.ctht.info/Table%2013%20DRB1%20DQA1%20DQB1%20associations%20in%20various%20populations.pdf
        if ($type eq "HLA-DQB1*05:01" && $$ref_dqa1[0] eq "HLA-DQA1*01:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*02:02" && $$ref_dqa1[0] eq "HLA-DQA1*02:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*03:04" && $$ref_dqa1[0] eq "HLA-DQA1*03:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*05:03" && $$ref_dqa1[0] eq "HLA-DQA1*01:04") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*06:01" && $$ref_dqa1[0] eq "HLA-DQA1*01:03") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*06:03" && $$ref_dqa1[0] eq "HLA-DQA1*01:03") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*06:04" && $$ref_dqa1[0] eq "HLA-DQA1*01:02") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*06:09" && $$ref_dqa1[0] eq "HLA-DQA1*01:02") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*06:02" && $$ref_dqa1[0] eq "HLA-DQA1*01:02") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*03:02" && $$ref_dqa1[0] eq "HLA-DQA1*03:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*04:02" && $$ref_dqa1[0] eq "HLA-DQA1*04:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*02:01" && $$ref_dqa1[0] eq "HLA-DQA1*05:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
        if ($type eq "HLA-DQB1*03:01" && $$ref_dqa1[0] eq "HLA-DQA1*05:01") {return $$ref_dqa1[0]."/".$type;} else {return $$ref_dqa1[1]."/".$type;}
      }
    }
  }else # DQA1 alleles cannot be typed, can only guess
  {
    if ($type eq "HLA-DQB1*05:01") {return "HLA-DQA1*01:01/".$type;}
    if ($type eq "HLA-DQB1*02:02") {return "HLA-DQA1*02:01/".$type;}
    if ($type eq "HLA-DQB1*03:04") {return "HLA-DQA1*03:01/".$type;}
    if ($type eq "HLA-DQB1*05:03") {return "HLA-DQA1*01:04/".$type;}
    if ($type eq "HLA-DQB1*06:01") {return "HLA-DQA1*01:03/".$type;}
    if ($type eq "HLA-DQB1*06:03") {return "HLA-DQA1*01:03/".$type;}
    if ($type eq "HLA-DQB1*06:04") {return "HLA-DQA1*01:02/".$type;}
    if ($type eq "HLA-DQB1*06:09") {return "HLA-DQA1*01:02/".$type;}
    if ($type eq "HLA-DQB1*06:02") {return "HLA-DQA1*01:02/".$type;}
    if ($type eq "HLA-DQB1*03:02") {return "HLA-DQA1*03:01/".$type;}
    if ($type eq "HLA-DQB1*04:02") {return "HLA-DQA1*04:01/".$type;}
    if ($type eq "HLA-DQB1*02:01") {return "HLA-DQA1*05:01/".$type;}
    if ($type eq "HLA-DQB1*03:01") {return "HLA-DQA1*05:01/".$type;}
  }

  return $type; 
}

