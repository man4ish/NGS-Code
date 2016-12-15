#!/usr/bin/perl
use Tie::IxHash;
tie (%hash, Tie::IxHash);
open(fh,"nucleotide.txt") or die "could not open the file\n";
@in=<fh>;
sub formatseq
{
    $len = length($_[0]);
    $qo = int ($len/70);
    $rem = $len%70;
    $i3=0;
    while($i3<$qo)
    {
         $pos = $i3*70;
         $le = 70;
         print  substr($_[0],$pos,$le),"\n"; 
         $i3++;
    }
   $lastpos = $len - $rem;
   print  substr($_[0],$lastpos,$rem),"\n";
}


for($i=0;$i<scalar(@in);$i++)
{
    if($in[$i]=~/>CCDS\d+\.\d+/)
    {
       $var= $&;
       chomp($var);
       push(@cdsid,$var);
       push(@cds,$i);
    }
}

for ($k=0;$k<scalar(@cds)-1;$k++)
{
     getseq($cds[$k],$cds[$k+1],$cdsid[$k]);
}


sub getseq
{

   $key  =  $_[2];
   $stpos=  $_[0];
   $endpos=  $_[1];
   $line = "";
   for($i1=$_[0]+1; $i1<$_[1] ;$i1++)
   {
       chomp($in[$i1]);
       $line .= $in[$i1];
   }
  $hash{$key} = $line;

open(fh,'CCDS.20080609.txt') or die "could not open the file\n";
while($input=<fh>)
{
      if ($input=~/(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t.*\t(\S)\t(\d+)\t(\d+)\t\[(.*)\]/)
      {
         @lendata="";
         @s=split(',',$8);   
         foreach $x(@s)
         {
                ($start,$end)=split('-',$x);
                 $leng = $end-$start;
                 push(@lendata,$leng);
         }
        $y = ">".$4;
      if( $hash{$y} ne '')
      {
          getsequence($y,\@lendata,\@s,$1,$2,$3,$5);
      }
    }
}

sub getsequence
{
    my $id = $_[0];
    my $sd=$_[2];
    my $chromosomeno = $_[3];
    my $geneaccession = $_[4];
    my $geneid = $_[5];
    my $chr_orient =$_[6];
    @se=@$sd;
    my $lea=$_[1];
    @ab = @$lea;
    $size = scalar(@ab);
    shift(@ab); 
   open(fh1,'info.txt') or die "could not oopen the file\n";
   %store =();
   while($input1=<fh1>)
   {
         if($input1=~/GeneID:(\d+)\tCDS\t(.*)\t(\w{2}\_\d+\.\d+)\t/)
         {
            $store{$1}=$3;
         }
   }
  my $sequence =$hash{$id};  
  for($ii=0;$ii<scalar(@ab);$ii++)
  {
      $ecs=$se[$ii];
      ($start,$end)=split('-',$ecs);
      $len=$ab[$ii]+1;
      $seq = substr($sequence,0,$len);
      $id=~s/>//;
      print ">GeneID:$geneid|seq|$geneaccession|Homo Sapiens Chromosome= $chromosomeno,chr_orient= $chr_orient,cds_start=   $start,cds_stop=$end,length=$len,transcript $store{$geneid},$id\n";
     formatseq($seq);
     $sequence=~s/$seq//;
  }
}






























