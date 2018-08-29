###########  detect neoantigens from a full length amino acid sequence  ############
# 
# prerequisite in path: 
# Rscript, iedb (MHC_I, MHC_II), python (python 2)
#
# input format:
# mhc_i, mhc_ii: folders to the iedb mhc1 and mhc2 binding prediction algorithms, http://www.iedb.org/
# fasta: amino acid sequence in fasta format
# type: HLA subtype, in the format of "allele class_allele 1_allele 2", e.g. "DRB1_HLA-DRB1*15:01_HLA-DRB1*12:01". Class A, B, C, DRB1, and DQB1 are allowed
# percentile_cutoff: percentile cutoff for binding affinity (0-100), recommended: 2
# output: output folder, safer to make "output" a folder that only holds results of this analysis job
#
#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($mhc_i,$mhc_ii,$fasta,$type,$percentile_cutoff,$output)=@ARGV;
my ($line,$path,$class,$len);
my %neoantigen_len=("I"=>[8..11],"II"=>[15]);

# prepare
$path=abs_path($0);
$path=~s/detect_neoantigen_aa\.pl//;
$path.="/script/";
unless (-d $output) {mkdir($output);}

# fake intermediate amino acid sequence files
foreach $class (keys %neoantigen_len)
{
  foreach $len (@{$neoantigen_len{$class}})
  {
    open(FILE_IN,$fasta);
    open(FILE_OUT,">".$output."/neoantigen_".$class."_".$len.".fa");

    $line=<FILE_IN>;
    print FILE_OUT ">protein-0-0-0-0-0-0\n";
    while ($line=<FILE_IN>) {print FILE_OUT $line;}

    close(FILE_IN);
    close(FILE_OUT);
  }
}

# fake typing file
open(FILE_OUT,">".$output."/typing.txt");
$type=~s/_/\t/g;
print FILE_OUT $type."\t0\n";
close(FILE_OUT);

# run affinity.pl
system_call("perl ".$path."/affinity.pl ".$mhc_i." ".$mhc_ii." ".$output."/neoantigen ".$output."/typing.txt ".$percentile_cutoff);
system_call("Rscript ".$path."/expressed_neoantigen_aa.R ".$output."/neoantigen ".$output." ".$path);

sub system_call
{
  my $command=$_[0];
  print "\n".$command."\n";
  system($command);
}

#perl /home2/twang6/software/immune/neoantigen/detect_neoantigen_aa.pl \
#/home2/twang6/software/immune/mhc_i \
#/home2/twang6/software/immune/mhc_ii \
#~/iproject/test/example.fa \
#"B_HLA-B*44:02_HLA-B*44:03" 2 \
#~/iproject/test/amino_acid

