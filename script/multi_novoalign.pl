#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use IO::Pipe;
use IO::Select;

my ($path,$thread,$index,$fastq1,$fastq2)=@ARGV;
my ($select,@select,$process,$pid,%pipes,$pm);
my ($first,@ready,$count,$line1,$line2);
my $command="-t 10 -o SAM -r all 100 -e 100 -i PE 200 140";
$SIG{CHLD}='IGNORE';

####  prepare for forking  ##############

$pm=Parallel::ForkManager->new($thread); # prepare for fork
$select=IO::Select->new(); # communication back to parent
foreach (1..$thread) {$pipes{$_}=IO::Pipe->new();} 

#------  child process starts  -------------------#

foreach $process (1..$thread) 
{	
  # do the fork
  $pid=$pm->start and next;
  $pipes{$process}->writer();
  $pipes{$process}->autoflush(1);

  # split fastq file
  open(FILE_IN1,$fastq1) or die "Cannot open fastq file to read!\n";
  open(FILE_IN2,$fastq2) or die "Cannot open fastq file to read!\n";
  open(FILE_OUT1,">".$fastq1."_".$process) or die "Cannot open fastq file to write!\n";
  open(FILE_OUT2,">".$fastq2."_".$process) or die "Cannot open fastq file to write!\n";
  $count=0;

  while ($line1=<FILE_IN1>)
  {
    foreach ((0..2)) {$line1.=<FILE_IN1>;}
    $line2=<FILE_IN2>;
    foreach ((0..2)) {$line2.=<FILE_IN2>;}

    if ($count++ % $thread==$process-1)
    {
      print FILE_OUT1 $line1;
      print FILE_OUT2 $line2;
    }
  }

  # novoalign
  system("novoalign -r None -t 10 -o Softclip -d ".$index." -f ".$fastq1."_".$process." ".$fastq2."_".
    $process." ".$command." > ".$path."/novo.sam_".$process);
  unlink($fastq1."_".$process);
  unlink($fastq2."_".$process);
  print {$pipes{$process}} $process;  
  print "Novoalign (".$process.") is done\n";
  $pm->finish; 
}

#---------  child process ends  ------------------#

#-------  parent process starts  ------------------#

# add child processes
foreach $process (1..$thread) 
{
  $pipes{$process}->reader();
  $select->add($pipes{$process});
}

# wait for all child processes
$first=1;
while ($select->count()>0)  
{
  @ready=$select->can_read(1);

  foreach (@ready)
  {
    $process=readline($_);
    $select->remove($_);
    if ($first==0)
    {
      system("grep -v -P \"^@\" ".$path."/novo.sam_".$process." | grep HLA > ".$path."/tmp.sam");
      system("cat ".$path."/novo.sam ".$path."/tmp.sam >".$path."/tmp2.sam");
      system("mv ".$path."/tmp2.sam ".$path."/novo.sam");
      unlink($path."/novo.sam_".$process);
      unlink($path."/tmp.sam");
      unlink($path."/tmp2.sam");
    }else
    {
      $first=0;
      system("mv ".$path."/novo.sam_".$process." ".$path."/novo.sam");
    }
  }
}

system("samtools view ".$path."/novo.sam -bS -h -F 4 > ".$path."/mpi.bam");
unlink($path."/novo.sam");
system("samtools sort -n -T ".$path."/tmp -O BAM -o ".$path."/mpi.sort.bam ".$path."/mpi.bam");

#-------  parent process ends  --------------------#

