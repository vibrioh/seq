#!/usr/bin/env perl
use Getopt::Long;

use DBI;
$PANDA=$ENV{'PANDA'};
$sample="";
$org="";
$method="";
$seqfile1="";
$seqfile2="";
#usage() if ( ! (GetOptions('help|?' => \$help, 'sample=s' => \$sample, 'org=s' => \$org, 'exp=s' => \$method, 'seq1=s' => \$seqfile1, 'seq2:s' => \$seqfile2)&& @ARGV<1)  or defined $help );

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



@refs=();
@cases=();
@samples=();
@experiments=();

if(-e "all_finshed.txt"){}else{
open(FIN, "sample.txt")||die;
if(defined($line=<FIN>)){}
$unfinished_count=0;
$eof=0;
while($eof==0){
  if($unfinished_count<5){
    if(defined($line=<FIN>)){
    chomp($line);
    @fields=split(/\t/, $line);
    $sample=$fields[1];
    $org=$fields[2];

    if($org eq "human"){

      $refGenome="$PANDA/db/bowtie/hg18.fa";
      $refAnnot="$PANDA/db/bowtie/hg18_RefSeq.gtf";
    }elsif($org eq "mouse"){

      $refGenome="$PANDA/db/bowtie/mm9.fa";
      $refAnnot="$PANDA/db/bowtie/mm9_RefSeq.gtf";
    }elsif($org eq "rat"){

      $refGenome="$PANDA/db/bowtie/rn4.fa";
      $refAnnot="$PANDA/db/bowtie/rn4_RefSeq.gtf";
    }


    $type=$fields[3];
    if($type eq "Reference"){
     push(@refs, $sample);

    }elsif($type eq "Case"){
     push(@cases, $sample);
    }
    $method=$fields[5];
    $seqfile1=$fields[6];
    $seqfile2=$fields[7];
    push(@samples, $sample);
    if(-e $sample."_tophat_finished.txt"){}else{
      system("perl $PANDA/code/tophat_cufflinks_single.pl --sample=".$sample." --exp=".$method." --org=".$org." --seq1=".$seqfile1."  --seq2=".$seqfile2." &");
      }
    }else{
      $eof=1;
   }
}
$finished_count=0;
foreach $sample (@samples){
  if(-e $sample."_tophat_finished.txt"){
    $finished_count++;
  }
}
$unfinished_count=($#samples+1)-$finished_count;
print "the numer of samples being  processed: ".$unfinished_count."\n";

sleep(60);

}
$finished_count=0;
foreach $sample (@samples){
  if(-e $sample."_tophat_finished.txt"){
    $finished_count++;
  }
}

while($finished_count<($#samples+1)){
sleep(1800);
$finished_count=0;
foreach $sample (@samples){
  if(-e $sample."_tophat_finished.txt"){
    $finished_count++;
  }
}
}


##########
#run cuffdiff for known transcripts
##########
$refStr="";
$i=0;
for $i (0 .. $#refs){
  if($refStr eq ""){
  $refStr=$refs[$i]."/accepted_hits.bam";
  }else{
  $refStr=$refStr.",".$refs[$i]."/accepted_hits.bam";
  }
}

#print $refStr."\n";
$i=0;
$caseStr="";
for $i (0 .. $#cases){
  if($caseStr eq ""){
  $caseStr=$cases[$i]."/accepted_hits.bam";
  }else{
  $caseStr=$caseStr.",".$cases[$i]."/accepted_hits.bam";
  }
}

#print $caseStr."\n";

#print "cuffdiff -p 4 -o known_Transcript_diff ".$refAnnot." ".$refStr." ".$caseStr."\n";
system("cuffdiff -p 5 -o known_Transcript_diff ".$refAnnot." ".$refStr." ".$caseStr); 


##########
#run cuffdiff for known transcripts
##########

 open(FOUT, ">assemblies.txt")||die;
 
$i=0;
for $i (0 .. $#refs){
   print FOUT $refs[$i]."/transcripts.gtf\n"; 
 }

$i=0;
 for $i (0 .. $#cases){
   print FOUT $cases[$i]."/transcripts.gtf\n";
 }
 close(FOUT);

 system("cuffmerge -p 5 -s ".$refGenome. " -g ".$refAnnot." assemblies.txt");

# system("cuffcompare -s ".$refGenome."  -r ".$refAnnot." merged_asm/merged.gtf");

 system("cuffdiff -p 5 -o novel_Transcript_diff merged_asm/merged.gtf ".$refStr." ".$caseStr);

 #system("R --no-save < $PANDA/code/combine_tophat.R");

 open(FOUT, ">all_finished.txt")||die;
 }

system("R --no-save < $PANDA/code/combine_data1.R");
