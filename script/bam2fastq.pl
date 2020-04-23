# convert one paired bam file to two fastq files
# prerequisite: samtools (>=1.4), sambamba
# bam_files: one bam file including path
# output: the output folder. Two fastq files will be written here: fastq1.fastq and fastq2.fastq. You can further gzip them
# thread: number of threads to use
#!/usr/bin/perl
use strict;
use warnings;

my ($bam_file,$output,$thread)=@ARGV;
my ($line,$name1,$name2,$read1,$read2,%reads);
my $count=0; # number of lines that have been written

# clean up
print "bam2fq for ".$bam_file."\n";
system("rm -f -r ".$output."/sambamba_tmp");
system("mkdir ".$output."/sambamba_tmp");
unlink($output."/sorted.bam");
  
# sambamba sort
print "First sambamba sort\n";
if (! -e $bam_file) 
{
  print $bam_file." doesn't exist!\n";
  exit;
}
system("sambamba sort --memory-limit=64GB --tmpdir=".$output."/sambamba_tmp -o ".$output."/sorted.bam -N -t ".$thread.
  " -l 0 -u ".$bam_file);

if (!-e $output."/sorted.bam") # if sambamba sort fails, possibly due to problem in bam file
{
  print "Warning: Second sambamba sort\n";
  system("sambamba view -f bam -h -v -l 0 -t ".$thread." ".$bam_file." -o ".$output."/sambamba_tmp/regenerate.bam");
  system("sambamba sort --memory-limit=64GB --tmpdir=".$output."/sambamba_tmp -o ".$output."/sorted.bam -N -t ".$thread.
    " -l 0 -u ".$output."/sambamba_tmp/regenerate.bam");
}

# bam2fq
print "Bam2fq\n";
system("samtools bam2fq -0 ".$output."/tmp0.fastq -1 ".$output."/tmp1.fastq -2 ".$output."/tmp2.fastq -n --threads ".$thread." ".$output."/sorted.bam");
unlink($output."/sorted.bam");
unlink($output."/tmp0.fastq");

print "Delete singletons\n";
open(FILE_IN1,$output."/tmp1.fastq");
open(FILE_IN2,$output."/tmp2.fastq");

open(FILE_OUT1,">".$output."/fastq1.fastq");
open(FILE_OUT2,">".$output."/fastq2.fastq");

while ($name1=<FILE_IN1>)
{
  $read1=$name1;
  foreach (0..2) {$line=<FILE_IN1>;$read1.=$line;}
  if (exists $reads{$name1})
  {
    print FILE_OUT1 $read1;
    print FILE_OUT2 $reads{$name1}; 
    $count++;
    delete $reads{$name1};
  }else
  {
    $reads{$name1}=$read1;
  }

  if (defined($name2=<FILE_IN2>))
  {
    $read2=$name2;
    foreach (0..2) {$line=<FILE_IN2>;$read2.=$line;}
    if (exists $reads{$name2})
    {
      print FILE_OUT2 $read2;
      print FILE_OUT1 $reads{$name2};
      $count++;
      delete $reads{$name2};
    }else
    {
      $reads{$name2}=$read2;
    } 
  }
}

while ($name2=<FILE_IN2>)
{
  $read2=$name2;
  foreach (0..2) {$line=<FILE_IN2>;$read2.=$line;}
  if (exists $reads{$name2})
  {
    print FILE_OUT2 $read2;
    print FILE_OUT1 $reads{$name2};
    $count++;
    delete $reads{$name2};
  }else
  {
    $reads{$name2}=$read2;
  }
}

close(FILE_OUT1);
close(FILE_OUT2);

close(FILE_IN1);
close(FILE_IN2);

# error checking
if ($count<100)
{
  print "Error: too few reads written!\n";
  unlink($output."/fastq1.fastq");
  unlink($output."/fastq2.fastq");
}

# clean up
unlink($output."/tmp1.fastq");
unlink($output."/tmp2.fastq");
END {system("rm -f -r ".$output."/sambamba_tmp");}
