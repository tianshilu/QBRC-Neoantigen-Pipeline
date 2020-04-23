#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/floor/;

# the first parameter is the index file
# the second parameter is the genome file
# the third parameter is the coordinates file
#   chr,strand,start,end,name (separated by tabs)
# the fourth parameter is the output sequence file
# the fifth parameter specifies the format of output
# 0: raw, 1: fasta, 2: stdout

my ($index,$genome,$coords,$seq,$format)=@ARGV;
my ($name,%index,@items,$chr,$strand,$start,$end,$s,$str,$line_s,$line_e,$sequence,$region);
my $jump=5000;

############  read index file ###############

open(FILE_IN,$index);

while (<FILE_IN>)
{
  $_=~s/\n//;
  @items=split("\t",$_);
  $index{$items[0]}=[$items[1],$items[2]];
}

close(FILE_IN);

############  query sequence  ##############

open(FILE_IN1,$coords);
open(FILE_IN2,$genome);
open(FILE_OUT,"> ".$seq);

while ($region=<FILE_IN1>)
{
  # locate start and end

  ($chr,$strand,$start,$end,$name)=split("\t",$region);
  $s=floor($start/$jump)*$jump+1; # this is the anchor
  if ($s==$start+1) {$s-=$jump;} # boundary case
  unless (exists $index{$chr."_".$s}) {print "The following region does not exist ".$region;}
  $sequence="";

  # the first segment

  @items=@{$index{$chr."_".$s}};
  seek FILE_IN2,$items[1],0; # go to the line
  $line_s=$s-$items[0]+1; # the coordinate of the first base on that line

  while ($line_s<=$start) # until the line in which $start resides is read
  {
    $str=<FILE_IN2>;
    unless (defined $str && $str!~/>/) {print "The following region does not exist ".$region;}
    $line_s+=(length($str)-1); # the coordinate of the first base on each line
  }

  $str=~s/\n//;
  $line_e=$line_s-1;
  $line_s-=length($str);

  if ($end>$line_e) # $start and $end are not on the same line
  {
    $sequence=substr $str,($start-$line_s);
  }else # $start and $end are on the same line
  {
    $sequence=substr $str,($start-$line_s),($end-$start+1);
    print_seq($sequence,$strand,$end-$start+1,$region,$format,$name);
    next;
  }

  # the middle segments

  $line_s+=length($str);

  while ($line_s<=$end)
  {
    $str=<FILE_IN2>;
    unless (defined $str && $str!~/>/) {print "The following region does not exist ".$region;}
    $str=~s/\n//;
    $line_s+=length($str);
    $sequence.=$str;
  }
  
  # the last segment
  
  if ($line_s-$end-1!=0) {$sequence=substr $sequence,0,-($line_s-$end-1);}
  print_seq($sequence,$strand,$end-$start+1,$region,$format,$name);
}

close(FILE_IN1);
close(FILE_IN2);
close(FILE_OUT);

exit;

sub print_seq
{
  my ($sequence,$strand,$len,$region,$format,$name)=@_;

  if ($len!=length($sequence)) {print "The following region does not exist ".$region;}

  if ($strand eq "-")
  {
    $sequence=~tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    $sequence=reverse($sequence);
  }

  $sequence=uc $sequence;
 
  if ($format==2)
  {
    print $sequence;
  }else
  {
    if ($format==1) {print FILE_OUT ">".$name;}
    print FILE_OUT $sequence."\n";
  }
}














