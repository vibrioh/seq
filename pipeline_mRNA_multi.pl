#!/usr/bin/env perl
use Getopt::Long;
use DBI;
$PANDA=$ENV{'PANDA'};
$sample="";
$org="";
$method="";
$seqfile1="";
$seqfile2="";
$splice_len=65;
$unique="Y";
$win=1000;
#usage() if ( ! (GetOptions('help|?' => \$help, 'sample=s' => \$sample, 'org=s' => \$org, 'exp=s' => \$method, 'seq1=s' => \$seqfile1, 'seq2:s' => \$seqfile2,'splice_len=i' => \$splice_len, 'win=i' => \$win, 'unique_align=s' => \$unique)&& @ARGV<1)  or defined $help );

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "usage: pipeline_mRNA.pl
              [--sample sample name]
              [--org species ( human, mouse, rat)]
              [--exp single end (SE) or paired-end (PE)
              [--seq1 read sequence 1 in single or paired end]
              [--seq2 read sequence 2 in paired-end (optional)]
              [--splice_len splice junction length at both ends [65]]
              [--win  window size to create density]
              [--help|? help]\n";
  exit;
}

@samples=();
@experiments=();
open(FIN, "sample.txt")||die;
if(defined($line=<FIN>)){}
$unfinished_count=0;
$eof=0;
while($eof==0){ 
  if($unfinished_count<6){
    if(defined($line=<FIN>)){
    chomp($line);
    @fields=split(/\t/, $line);
    $sample=$fields[1];
    $org=$fields[2];
    $method=$fields[5];
    $seqfile1=$fields[6];
    $seqfile2=$fields[7];
    $splice_len=int $fields[8];
    $win=int $fileds[9];
    $unique= $fields[10];
    push(@samples, $sample);
    if(-e $sample."_finished.txt"){}else{
    system("perl $PANDA/code/pipeline_mRNA_single.pl --sample=".$sample." --exp=".$method." --org=".$org." --seq1=".$seqfile1."  --seq2=".$seqfile2."  --win=".$win." --splice_len=".$splice_len." --unique_align=".$unique." &");
    }
  }else{
    $eof=1;
   }
}

$finished_count=0;
foreach $sample (@samples){
  if(-e $sample."_finished.txt"){
    $finished_count++;
  }
}
$unfinished_count=($#samples+1)-$finished_count;
print "the numer of samples being  processed: ".$unfinished_count."\n";

sleep(60);

}
$finished_count=0;
foreach $sample (@samples){
  if(-e $sample."_finished.txt"){
    $finished_count++;
  }
}

while($finished_count<($#samples+1)){
sleep(1800);
$finished_count=0;
foreach $sample (@samples){
  if(-e $sample."_finished.txt"){
    $finished_count++;
  }
}
}
system("R --no-save < $PANDA/code/combine_data1.R");
          
