#!/usr/bin/perl
use strict;
use warnings;

# the first parameter is the genome fasta file
# the second parameter is the output index file

my ($chr,$i,@items,$len,$pos,$next,$pre_len);
my $jump=5000;

open(FILE_IN,$ARGV[0]);
$pos=tell FILE_IN;
open(FILE_OUT,">".$ARGV[1]);

while (<FILE_IN>)
{
  if ($_=~/>(.+)\n/) 
  {
    @items=split("\t",$1);
    $chr=$items[0];
    $len=0;
    $next=0;
  }else
  {
    $pre_len=$len;
    $len+=(length($_)-1);

    if ($len-$jump*$next>0)    
    {
      print FILE_OUT $chr."_".($jump*$next+1)."\t".($jump*$next+1-$pre_len)."\t".$pos."\n";
      $next++;
    }
  }

  $pos=tell FILE_IN;
}

close(FILE_IN);
close(FILE_OUT);

exit;
