# prerequisite in path: novoalign, samtools, Athlates (need lib64 of gcc>=5.4.0 in LD_LIBRARY_PATH)
# fastq1, fastq2: fastq files for exome-seq
# output_dir: output directory
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use IO::Pipe;
use IO::Select;

my ($fastq1,$fastq2,$output_dir,$script_path,$thread)=@ARGV;
my ($key,$line,$bed,$allele,$bam1,$bam2,$typing_folder,@items,%inferred);
my ($select,$pipe,%pipes,$pm,$pid);
my $count=0;

# find home directory of typing
open(TP, "which typing |") or die "Error: Athlates is not in home path!";
$typing_folder=<TP>;
$typing_folder=~s/bin\/typing\n//;
close(TP);

unless (-e $fastq1) {die("Error: ".$fastq1." doesn't exist\n");}
unless (-e $fastq2) {die("Error: ".$fastq2." doesn't exist\n");}
system_call("perl ".$script_path."/multi_novoalign.pl ".$output_dir." ".$thread." ".$typing_folder."/db/ref/ref.nix ".
  $fastq1." ".$fastq2);

# process each allele
$SIG{CHLD}='IGNORE';
$pm=Parallel::ForkManager->new(9);
$select=IO::Select->new(); # communication back to parent

foreach $allele (("DRB1","DRB3","DRB4","DRB5","DQA1","DQB1","A","B","C"))
{
  # initiate child process
  $pipes{$allele}=IO::Pipe->new();
  $pid = $pm -> start and next;
  $pipes{$allele}->writer();
  $pipes{$allele}->autoflush(1);

  system_call("rm -f -r ".$output_dir."/".$allele);
  system_call("mkdir ".$output_dir."/".$allele);

  # extract reads in or outside of HLA alleles
  foreach $bed (("hla.".$allele.".bed","hla.non-".$allele.".bed"))
  {
    system_call("samtools view -b -L ".$typing_folder."/db/bed/".$bed." ".$output_dir."/mpi.sort.bam > ".$output_dir."/".$allele."/tmp.bam");
    system_call("samtools view -h -o ".$output_dir."/".$allele."/tmp.sam ".$output_dir."/".$allele."/tmp.bam");
    system_call("grep -v -P \"^@\" ".$output_dir."/".$allele."/tmp.sam >".$output_dir."/".$allele."/tmp_body.sam");
    system_call("grep -P \"^@\" ".$output_dir."/".$allele."/tmp.sam >".$output_dir."/".$allele."/tmp_header.sam");
    system_call("sort -k 1,1 -k 3,3 ".$output_dir."/".$allele."/tmp_body.sam >".$output_dir."/".$allele."/tmp_body_sort.sam");
    system_call("cat ".$output_dir."/".$allele."/tmp_header.sam ".$output_dir."/".$allele."/tmp_body_sort.sam >".$output_dir."/".$allele."/tmp_sort.sam");
    system_call("samtools view -bS ".$output_dir."/".$allele."/tmp_sort.sam >".$output_dir."/".$allele."/".$bed.".bam");
  }

  # typing
  $bam1=$output_dir."/".$allele."/hla.".$allele.".bed.bam";
  $bam2=$output_dir."/".$allele."/hla.non-".$allele.".bed.bam";
  system_call($typing_folder."/bin/typing -bam ".$bam1." -exlbam ".$bam2." -hd 5 -msa ".
    $typing_folder."/db/msa/".$allele."_nuc.txt -o ".$output_dir."/".$allele."/typing");

  if (!-e $output_dir."/".$allele."/typing.typing.txt")
  {
    print "Warning: Allele ".$allele." cannot be typed(1)!\n";
    print {$pipes{$allele}} "NT";
    $pm->finish;
  }

  # skip unuseful lines in report
  open(FILE_IN,$output_dir."/".$allele."/typing.typing.txt"); 
  while ($line=<FILE_IN>)
    {if ($line=~/Inferred Allelic Pairs/) {last;}}
  $line=<FILE_IN>;

  # extract typed results
  %inferred=();
  while ($line=<FILE_IN>)
  {
    $line=~s/\n//;
    @items=split("\t\t",$line);
    $items[0]=~/(.*\*[0-9]+\:[0-9]+)/;
    $items[0]="HLA-".$1;
    $items[1]=~/(.*\*[0-9]+\:[0-9]+)/;
    $items[1]="HLA-".$1;
    $inferred{$allele."\t".$items[0]."\t".$items[1]}=$items[2];  
  }
  close(FILE_IN);

  # send output to parent
  if (scalar(keys %inferred)>=1) 
  {
    $key=(keys %inferred)[0];
    print "Typed: ".$key."\t".$inferred{$key}."\n";
    print {$pipes{$allele}} $key."\t".$inferred{$key}."\n";
  }else
  {
    print "Warning: Allele ".$allele." cannot be typed(2)!\n";
    print {$pipes{$allele}} "NT";
  }
  $pm->finish;
}

# set up communications
foreach $allele (keys %pipes) 
{
  $pipes{$allele}->reader();
  $select->add($pipes{$allele});
}

# get output from children
open(FILE_OUT,">".$output_dir."/typing.txt") or die "Error: Cannot write to typing file!\n";

while ($select->count()>0)
{
  foreach ($select->can_read(1))
  {
    $line=readline($_);
    $select->remove($_);
    if ($line ne "NT") 
    {
      $count++;
      print FILE_OUT $line;
    }
  }
}

close(FILE_OUT);

# clean up
foreach $allele (keys %pipes) {system_call("rm -f -r ".$output_dir."/".$allele);}
unlink($output_dir."/mpi.bam");
unlink($output_dir."/mpi.sort.bam");

unlink($output_dir."/error_no_hla_allele_can_be_typed.txt");
if ($count==0) 
{
  unlink($output_dir."/typing.txt");
  system("echo \"\" > ".$output_dir."/error_no_hla_allele_can_be_typed.txt");
  die("Warning: No HLA allele can be typed! Halt execution!\n");
}

sub system_call
{
  my $command=$_[0];
  print $command."\n";
  system($command);
}

