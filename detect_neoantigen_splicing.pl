# detect neoantigens from somatic splicing events
# output from this pipeline still needs to be examined manually on a case by case basis
# may need to set "ulimit -n 50000" in .bash_profile
# need >=256GB of memory
#
# prerequisite in path: novoalign, Athlates (need lib64 of gcc>=5.4.0 in LD_LIBRARY_PATH), samtools (>=1.4), sambamba
# stringtie (>=1.3.2), STAR (>=2.7.2), blast+, python (python 2)
# 
# build: genome build version, hg38 or hg19
# reference: path to reference genome 
# fastq1,fastq2,fastq3,fastq4: fastq files of the normal and tumor RNA-seq data for HLA typing and stringtie analysis. 
#        Or bam files. refer to detect_neoantigen.pl for format
# star: STAR index. already generated as part of the reference genome data bundle
# output: output folder
# thread: number of thread
# mhc_i, mhc_ii: folders to the iedb mhc1 and mhc2 binding prediction algorithms, http://www.iedb.org/
# percentile_cutoff: percentile cutoff for binding affinity (0-100), recommended: 2
# hla_file: optional. a pre-existing HLA typing file. default is NA. if specified, will skip the HLA typing process
#
#perl /home2/twang6/software/immune/neoantigen/detect_neoantigen_splicing.pl \
#hg38 /home2/twang6/data/genomes/hg38/hs38d1.fa \
#/project/bioinformatics/Xiao_lab/shared/genomics/IL2/RNA_seq_latest/SHI181-27-8134_N_RNA_wholernaseq-24-19119_S7_R1_001.fastq.gz \
#/project/bioinformatics/Xiao_lab/shared/genomics/IL2/RNA_seq_latest/SHI181-27-8134_N_RNA_wholernaseq-24-19119_S7_R2_001.fastq.gz \
#/project/bioinformatics/Xiao_lab/shared/genomics/IL2/RNA_seq_latest/SHI181-27-8122_T_RNA_wholernaseq-24-19119_S6_R1_001.fastq.gz \
#/project/bioinformatics/Xiao_lab/shared/genomics/IL2/RNA_seq_latest/SHI181-27-8122_T_RNA_wholernaseq-24-19119_S6_R2_001.fastq.gz \
#/home2/twang6/data/genomes/hg38/STAR \
#/project/SCCC/Wang_lab/shared/tmp 32 \
#/project/bioinformatics/Xiao_lab/shared/neoantigen/code/mhc_i \
#/project/bioinformatics/Xiao_lab/shared/neoantigen/code/mhc_ii \
#2 NA
#
#!/usr/bin/perl
use warnings;
use strict;
use Cwd 'abs_path';
use Parallel::ForkManager;

###########  prepare  #############################

my ($build,$reference,$fastq1,$fastq2,$fastq3,$fastq4,$star,$output,$thread,
  $mhc_i,$mhc_ii,$percentile_cutoff,$hla_file)=@ARGV;
my $debug=1;
my $fpkm_cutoff=2;

my ($type,$path,$bam_file,$line,@seqs,@items,$chr,$transcript_id,$fpkm,$strand);
my ($gene_id,$to_print,$reference_folder,$a,$b,%white_list,%final_protein);
my ($pid,$pm,$len,$class,$protein,$wt_proteins,$normal_proteins,$i,$j,$epitope);
my (%fastqs,$stringtie_parameters,$k);

my %neoantigen_len=("I"=>[8..11],"II"=>[15]); # length of neoantigen for HLA A/B/DRB/DQB
%fastqs=("normal"=>[$fastq1,$fastq2],"tumor"=>[$fastq3,$fastq4]);

$path=abs_path($0);
$path=~s/detect_neoantigen_splicing\.pl//;

$reference=~/(.*\/).*$/;
$reference_folder=$1;

system_call("rm -f -r ".$output."/binders");
unless (-d $output) {mkdir($output) or die "Error: cannot create directory ".$output."!\n";}
if (!-w $output) {die "Error: directory ".$output." is not writable!\n";}

##########  HLA typing and stringtie  ##########

$pm=Parallel::ForkManager->new(2);

foreach $type (keys %fastqs)
{
  $pid = $pm -> start and next;

  $fastq1=${$fastqs{$type}}[0];
  $fastq2=${$fastqs{$type}}[1];
  system_call("rm -f -r ".$output."/".$type);
  system_call("mkdir ".$output."/".$type);

  # HLA typing  
  if ($fastq1 eq "bam") # given bam files 
  {
    system_call("perl ".$path."/script/bam2fastq.pl ".$fastq2." ".$output."/".$type." ".$thread);
    if (! -e $output."/".$type."/fastq1.fastq")
    {
      print "Error: Fastq file doesn't exist!\n";
      exit;
    }
  }else  # given fastq files 
  {
    system_call("cp ".$fastq1." ".$output."/".$type."/fastq1.fastq.gz");
    system_call("gzip -d -f ".$output."/".$type."/fastq1.fastq.gz");

    system_call("cp ".$fastq2." ".$output."/".$type."/fastq2.fastq.gz");
    system_call("gzip -d -f ".$output."/".$type."/fastq2.fastq.gz");
  }

  if ($type eq "tumor")
  {
    if ($hla_file eq "NA")
    {
      system_call("perl ".$path."/script/hla.pl ".$output."/".$type."/fastq1.fastq ".
        $output."/".$type."/fastq2.fastq ".$output."/".$type." ".$path."/script ".$thread);
      system_call("cp ".$output."/".$type."/typing.txt ".$output);
    }else
    {
      system_call("cp ".$hla_file." ".$output."/typing.txt"); 
    }
  }

  #stringtie  
  system_call("STAR --runThreadN ".$thread." --genomeDir ".$star." --readFilesIn ".$output."/".$type."/fastq1.fastq ".
    $output."/".$type."/fastq2.fastq --outFileNamePrefix ".$output."/".$type."/ --outSAMtype BAM SortedByCoordinate ".
    "--outSAMattributes NH HI AS nM NM XS");
  unlink($output."/".$type."/fastq1.fastq");
  unlink($output."/".$type."/fastq2.fastq");
  $bam_file=$output."/".$type."/Aligned.sortedByCoord.out.bam";

  if ($type eq "tumor") {$stringtie_parameters="-f 0.15 -j 5 -c 10 -M 0.7";} else {$stringtie_parameters="-f 0.05 -j 2";}
  system_call("stringtie ".$bam_file." -o ".$output."/".$type."/out.gtf -G ".$reference_folder."/".
    $build."_genes.gtf -p ".$thread." -A ".$output."/".$type."/gene_abundance.txt ".$stringtie_parameters);
  unlink($bam_file);
 
  # obtain DNA sequences of novel transcripts 
  open(FILE,$output."/".$type."/out.gtf");
  open(COORD,">".$output."/".$type."/coordinates.txt");

  $line=<FILE>;
  $line=<FILE>;
  $transcript_id=$strand=$chr=$to_print="";
  $fpkm=0;

  while ($line=<FILE>)
  {
    if ($line=~/gene_id\s\"(.*?)\";.*reference_id\s\"(.*?)\";/) # known transcript/gene
    {
      $white_list{$1}=$2; # only keep transcripts from genes with at least one known transcript. This is alt splicing
      next;
    } 
  
    unless ($line=~/^chr.*\t[\+\-]\t/) {next;} # need to know strand and need to start with chr (regular chromosome)
    if ($line=~/chrM/) {next;} # the pipeline current does not work with chrM transcripts. may want to change this behavior
                               # stringtie sometimes gave out of boundary errors with chrM

    if ($line=~/(.*)\tStringTie\ttranscript.*\t([\+\-])\t.*gene_id\s\"(.*?)\";.*transcript_id\s\"(.*?)\";\scov.*FPKM\s\"(.*?)\";\sTPM/)
    {
      if ($4 ne $transcript_id)
      {
        if ($fpkm>$fpkm_cutoff || $type eq "normal") {print COORD $to_print;} # more lenient on normal samples
        $chr=$1;
        $strand=$2;
        $gene_id=$3;
        $transcript_id=$4;
        $fpkm=$5;
        $to_print="";
      }
    }else
    {
      $line=~/exon\t([0-9]+)\t([0-9]+)/;
      $to_print.=$chr."\t".$strand."\t".$1."\t".$2."\t".$transcript_id."\n";
    }
  }

  if ($fpkm>$fpkm_cutoff || $type eq "normal") {print COORD $to_print;}
  close(FILE);
  close(COORD);

  system_call("perl ".$path."/script/query.pl ".$reference_folder."/index.text ".$reference." ".
    $output."/".$type."/coordinates.txt ".$output."/".$type."/sequences.txt 0");

  # obtain protein sequences 
  $a=0;open(FILE,$output."/".$type."/sequences.txt");while (<FILE>) {$a++;};close(FILE);
  $b=0;open(FILE,$output."/".$type."/coordinates.txt");while (<FILE>) {$b++;};close(FILE);
  if ($a!=$b) {die $type.": Not all coordinates are valid!";}

  # get full DNA sequences
  open(COORD,$output."/".$type."/coordinates.txt");
  open(SEQ,$output."/".$type."/sequences.txt");
  open(PRO,">".$output."/".$type."/full_dna.txt");
  $transcript_id=$gene_id="";
  @seqs=();

  while ($line=<COORD>)
  {
    @items=split("\t",$line);

    if ($transcript_id ne $items[4])
    {
      print PRO full_dna(\@seqs,$transcript_id,$strand,$gene_id,\%white_list);
      @seqs=();
      $transcript_id=$items[4];
      $strand=$items[1];
      $gene_id=$transcript_id;
      $gene_id=~s/\.[0-9]*?\n//;
    }

    $line=<SEQ>;
    $line=~s/\n//;
    push @seqs,$line;
  }

  print PRO full_dna(\@seqs,$transcript_id,$strand,$gene_id,\%white_list);
  close(COORD);
  close(SEQ);
  close(PRO);

  # predict CDS and protein sequences
  system_call("rm -f -r ".$output."/".$type."/cds_search*");
  $ENV{'BLASTDB'}=$path."/data/blast/swissprot";
  system_call("cd ".$output."/".$type.";".$path."/script/TransDecoder/TransDecoder.LongOrfs -t ".$output."/".$type.
    "/full_dna.txt -m 50 -S --output_dir ".$output."/".$type."/cds_search");
  system_call("blastp -num_threads ".$thread." -query ".$output."/".$type."/cds_search/longest_orfs.pep -db ".$path.
    "/data/blast/swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > ".$output."/".$type."/cds_search/blastp.outfmt6");
  system_call("cd ".$output."/".$type.";".$path."/script/TransDecoder//TransDecoder.Predict -t ".$output."/".$type.
    "/full_dna.txt --single_best_only --retain_blastp_hits ".$output."/".$type."/cds_search/blastp.outfmt6 --output_dir ".
    $output."/".$type."/cds_search");

  # find good ORFs
  %final_protein=();
  open(ORF,$output."/".$type."/cds_search/longest_orfs.cds.best_candidates.gff3.revised_starts.gff3");

  while ($line=<ORF>)
  {
    unless ($line=~/\tCDS\t.*ID\=cds\.(.*p[0-9]+);Parent/) {next;}
    $final_protein{$1}=1;
  }

  close(ORF);

  # extract the good ones' protein sequences
  open(PEP,$output."/".$type."/cds_search/longest_orfs.pep"); # extract protein sequences of good ORFs

  while ($line=<PEP>)
  {
    $line=~/>(.*?)\stype/;
    $transcript_id=$1; # this a mis-namer. this should be protein id
    $line=<PEP>;
    $line=~s/\*\n//;
    $line=~s/\n//;
    unless (exists $final_protein{$transcript_id}) {next;}
    $final_protein{$transcript_id}=$line;
  }

  close(PEP);
 
  # write down the good protein sequences 
  open(FILE,">".$output."/".$type."/final_protein.txt");
  foreach (keys %final_protein) {print FILE $_." ".$final_protein{$_}."\n";}
  close(FILE);

  #finish up
  $pm->finish;
}

$pm->wait_all_children;

#############  get wt and normal proteins  #############

# get normal proteins
open(FILE,$output."/normal/final_protein.txt");
$normal_proteins="";

while ($line=<FILE>)
{
  $line=~s/\n//;
  @items=split(" ",$line);
  $normal_proteins.=($items[1]." ");
}

close(FILE);

# get tumor proteins
open(FILE,$output."/tumor/final_protein.txt");
%final_protein=();

while ($line=<FILE>)
{
  $line=~s/\n//;
  @items=split(" ",$line);
  $final_protein{$items[0]}=$items[1];
}

close(FILE);

# get wt proteins
open(PT,$path."/data/human_protein_seq.fasta") or die "Error: Cannot open the human protein sequence file!\n";

while ($line=<PT>)
{
  if ($line=~/>/)
  {
    $wt_proteins.=" ";
  }else
  {
    $line=~s/\n//;
    $wt_proteins.=$line;
  }
}

close(PT);

############  predict neoantigens  ######################

system_call("rm -f -r ".$output."/binders");
system_call("mkdir ".$output."/binders");
$pm=Parallel::ForkManager->new(5);

foreach $class (("I","II"))
{
  foreach $len (@{$neoantigen_len{$class}})
  {
    $pid = $pm -> start and next;
    print "Working on class ".$class." length ".$len."\n";
    open(FILE,">".$output."/binders/neoantigen_".$class."_".$len.".fa");

    foreach $transcript_id (keys %final_protein)
    {
      $i=0;
      $protein=$final_protein{$transcript_id};

      while ($i<=length($protein)-$len)
      {
        $j=look_ahead($i,$protein,\$wt_proteins,\$normal_proteins,$len); # can safely jump this far
        if ($j<0) # don't jump. this position will lead to a new epitope with length at least $len
        {
          # very important here. if the protein sequence is new from the beginning. 
          # this is extremely unlikely to be an alt splicing event
          # most likely an error in annotation and searching of orf
          if ($i<=5 || $i+$len-1>=length($protein)) {last;}

          $k=1; # extend the epitope as long as possible
          while ($i+$k+$len-1<=length($protein)-1)
          { 
            $epitope=substr $protein,$i+$k,$len;
            if (index($wt_proteins,$epitope)==-1 && 
                  index($normal_proteins,$epitope)==-1) {$k++;} else {last;}
          }

          $epitope=substr $protein,$i,$len+$k-1;
          print FILE ">".$transcript_id."-".$i."-".($i+$len-1+$k-1)."-NA-NA-NA-NA\n".$epitope."\n";
          $i+=($k+1);
        }else
        {
          $i+=$j+1;
        }   
      }
    }    

    close(FILE);
    $pm->finish;
  }
}

$pm->wait_all_children;

system_call("perl ".$path."/script/affinity.pl ".$mhc_i." ".$mhc_ii." ".$output."/binders/neoantigen ".
  $output."/typing.txt ".$percentile_cutoff);
system_call("Rscript ".$path."/script/splicing_cleanup.R ".$output."/binders/neoantigen ".$output." ".$path."/script ".$debug);

############  cleanup  ##################

system_call("cp ".$output."/tumor/coordinates.txt ".$output);
system_call("cp ".$output."/tumor/final_protein.txt ".$output);
system_call("cp ".$output."/tumor/full_dna.txt ".$output);
system_call("cp ".$output."/tumor/out.gtf ".$output."/tumor.gtf");
system_call("cp ".$output."/normal/out.gtf ".$output."/normal.gtf");
system_call("cp ".$output."/tumor/sequences.txt ".$output);
system_call("cp ".$output."/tumor/cds_search/longest_orfs.pep ".$output);
system_call("cp ".$output."/tumor/cds_search/longest_orfs.cds.best_candidates.gff3.revised_starts.gff3 ".$output);

if ($debug==0)
{
  system_call("rm -f -r ".$output."/normal");
  system_call("rm -f -r ".$output."/tumor");
  system_call("rm -f -r ".$output."/binders");
}

#############  subroutines  ############

# can move the index by this amount of distance with no problem, because the covered sequences are all wild type
sub look_ahead 
{
  my ($i,$protein_seq,$wt_proteins_ref,$normal_proteins_ref,$len)=@_;
  my ($temp,$j,$j0);
  my @look_ahead_dists=(0,10,20,50,200,500,1000,2000,4000,8000);

  $j0=$#look_ahead_dists;
  foreach $j (0..$#look_ahead_dists)
  {
    if ($i+$len+$look_ahead_dists[$j]-1>=length($protein_seq)) 
    {
      $temp=substr $protein_seq,$i;
      if (index($$normal_proteins_ref,$temp)!=-1 || 
            index($$wt_proteins_ref,$temp)!=-1) # all the rest is WT
      {
        return length($protein_seq); # just move to the end
      }else
      {
        $j0=$j-1;
        last;
      }
    }
    $temp=substr $protein_seq,$i,$len+$look_ahead_dists[$j];
    if (index($$normal_proteins_ref,$temp)==-1 &&
          index($$wt_proteins_ref,$temp) == -1) 
    {
      $j0=$j-1;
      last;
    }
  }

  if ($j0==-1) {return -1;} else {return $look_ahead_dists[$j0];}
}

sub full_dna
{
  my ($seqs_ref,$transcript_id,$strand,$gene_id,$white_list_ref)=@_;
  my (@seqs,$dna,$protein);
  if ($transcript_id eq "") {return "";}
  unless (exists $white_list_ref->{$gene_id}) {return "";}

  # concatenate DNA sequences
  @seqs=@$seqs_ref;
  if ($strand eq "-") {@seqs=reverse(@seqs);}
  $dna=join("",@seqs);

  # identify ORF
  ">".$transcript_id.$dna."\n"; 
}

sub system_call
{
  my $command=$_[0];
  print "\n".$command."\n";
  system($command);
}


