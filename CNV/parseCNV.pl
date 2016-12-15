#!/usr/bin/perl -w
open (fh , $ARGV[0]) or die "could not open file\n";
print "GenomeID\tUnique\tChromosome\tOffset\tVariant Type\tNucleotide_Variant - single base\tSupporting Reads (forward strand)\tSupporting Reads (reverse strand)\tSupporting Reads (inverse strand)\tValue: Total_Read_Density or CN\tAvg_QScore\tdbSNP ID\tGene_Name\tmiRNA ID\tPromoter\tExon\tZygosity\tLocal_Copy_Number\tKnown_CNVRegion\tTransversions(Tv)_Transitions(Ts)\tSynonymous(S)_Non-Synonymous(NS)\tProtein_Variant\tMissense\tNonsense\tSNV Confident Score\tSV - Length or CN - end of offset\tSV Link ID\tSV Edge Support\n";
while($input = <fh>)
{
@vec = split("\t",$input);
chomp($vec[8]);
$a = substr($vec[8],0,length($vec[8])-1);
$b = $vec[3];
$c= $vec[6]; 
$d = $vec[4]; 
$e = $vec[7]; 
$f = $vec[5]; 
$g = $vec[2];
  if($input !~ /#Chr/)
  {
  if($e eq "")
  {
    $e = "-";
  }
if($e eq "")
  {
    $e = "-";
  }
 
  if($b eq "")
  {
    $b = "-";
  }
  if($c eq "")
  {
    $c = "-";
  }
  if($d eq "")
  {
    $d = "-";
  }
  if($f eq "")
  {
    $f = "-";
  }
  if($e eq "")
  {
    $e = "-";
  }
if($g eq "")
  {
    $g = "-";
  }

  print  "$ARGV[1]\t$ARGV[2]\t$vec[0]\t$vec[1]\t$ARGV[3]\t-\t-\t-\t-\t$b\t-\t$c\t$d\t$e\t$a\t$f\t-\t-\t-\t-\t-\t-\t-\t-\t-\t$vec[2]\t-\t-\n";


}

}
