#!/usr/bin/perl 

#######################################################################################
#perl ftest.pl hs_ref_chr1.fa 1                                                       #
#######################################################################################
open(out,">>WholeExon_Database.txt") or die "PBLM!";
sub splitsequence {
   my $seque = $_[0];
   my $lendaata= $_[1];
   my $geneaa = $_[2];
   my $orientaa = $_[3];
   my $transcriptaa = $_[4];
   my $chrnoaa = $_[5];
   @a = @$lendaata;   
  for($iii=0;$iii<scalar(@a);$iii++) {
    ($asd1,$asd2)=split('-',$a[$iii]);
     $len=$asd2-$asd1+1;
     $seqio = substr($seque,0,$len-1);
     print out ">$geneaa|ref|$transcriptaa|, Homosapiens chromosome = $chrnoaa, cds start = $asd1, cds stop = $asd2, length = $len, chr_orient = $orientaa\n"; 
     #print out "$seqio\n";
      @sdfgh1 = split("",$seqio);
for($iiii1=1;$iiii1<=scalar(@sdfgh1);$iiii1++){
print out "$sdfgh1[$iiii1]";
if($iiii1%70==0) {print out "\n";}
}print out "\n";
     $seque=~s/$seqio//;
  }
}

sub rev {
@revcom="";
 @s= split('',$_[0]);
 @t =reverse(@s);
 for ($ik=0;$ik<scalar(@t);$ik++){
   if($t[$ik] eq 'A'){push(@revcom,'T');}
   if($t[$ik] eq 'C'){push(@revcom,'G');}
   if($t[$ik] eq 'G'){push(@revcom,'C');}
   if($t[$ik] eq 'T'){push(@revcom,'A');}
  }
$seqd = join('',@revcom);
return $seqd;
}


open(hs,$ARGV[0]) or die "could not open file";

 $y='';
 $i=0;
  while($seq=<hs>) {
     if($seq!~/^>gi/) {
      chomp($seq);
      $y .= $seq;
      push(@line,$seq);
     }
      $i++;
   }


sub getfisrtline {
  $line=$_[0];
  @base = split('',$line);
  $firstline ='';
   for ($c=$_[1]; $c < scalar(@base);$c++) {
     $firstline .= $base[$c];  
    }
return $firstline;
}

sub getlastline {
$line=$_[0];
  @base = split('',$line);
  $lastline ='';
   for ($c=0; $c < $_[1];$c++) {
    $lastline .= $base[$c];  
    }
return $lastline;
}

sub getbet {
@sdf=split('',$_[2]);
$sr='';
  for($k=$_[0];$k<$_[1];$k++) {
    $sr .= $sdf[$k];
   }
return $sr;
}

sub getseq {
  $start_line  =  $_[0]/70;
  $stop_line  =  $_[1]/70;
  $startline = int ($start_line );
  $remainder_start = $_[0]%70;
  $remainder_stop = $_[1]%70;
  $stopline = int ($stop_line);
  $output='';
   for($id = $startline+1; $id < $stopline ;$id++  )  #getting all the lines in between the start and stop line
    {
     $output .= $line[$id];    
    }
    
    if($startline eq $stopline)
      {
       $seqout= getbet($remainder_start,$remainder_stop,$line[$startline]);
      return $seqout;
      }
     elsif ($startline eq $stopline-1) {
         $firstseq=getfisrtline($line[$startline], $remainder_start);
         $lastseq =getlastline($line[$stopline], $remainder_stop);
         $seqout=$firstseq.$lastseq;
       return $seqout;
       }
else {
$firstseq=getfisrtline($line[$startline], $remainder_start);
$lastseq =getlastline($line[$stopline], $remainder_stop);
$seqout=$firstseq.$output.$lastseq;
return $seqout;
}
}

sub joinseq {
$genea = $_[0];
$westa = $_[1];
$orienta = $_[2];
$transcripta =$_[3];
$chrnoa = $_[4];
@west = @$westa;
$tyseq = "";
for($kk = 0;$kk < scalar(@west); $kk++) {
  ($mm,$nn) = split('-',$west[$kk]);
  $mmnnseq= getseq($mm,$nn);
  $tyseq .= $mmnnseq; 
 }
$revtyseq = rev($tyseq);
splitsequence($revtyseq,\@west,$genea,$orienta,$transcripta,$chrnoa); 
}
%orient = (); 
%transcript = ();
%chrno = ();
open(fk,"info.txt") or die "PBLM!";
while($input=<fk>) {
 if($input=~/\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t\S+\t(\S+)\t(\S+)/) {
   $transcript{$5} = $7;
   $y=$2."-".$3;
   push(@geneid,$5);
   $chrno{$5} = $1 ;
   $orient{$5} = $4; 
   $hsd{$5} = "";
   push(@offset,$y);
  }
}
%geneidvalue = ();
foreach $id(keys %hsd) {
   @value= ();
   for($i=0;$i<scalar(@geneid);$i++) {
        if($id eq $geneid[$i] ) {
         push(@value,$offset[$i]);      
         }
     }
   $geneidvalue{$id} = [@value] ;
 }

while (($kkk, $vvv) = each %geneidvalue) {
if($chrno{$kkk} == $ARGV[1])
{
if($orient{$kkk} eq '-') {
   joinseq($kkk,$vvv,$orient{$kkk},$transcript{$kkk},$chrno{$kkk});         
} 
 else {
    @abcd = @$vvv;
    for($dd=0; $dd< scalar(@abcd); $dd++) {
       ($m,$n) = split('-',$abcd[$dd]);
       $sseeqq = getseq($m,$n);
        $lleenn = $n - $m + 1;
       print out ">$kkk|ref|$transcript{$kkk}|, Homosapiens chromosome = $chrno{$kkk}, cds start = $m, cds stop = $n, length = $lleenn, chr_orient = $orient{$kkk}\n"; 
     @sdfgh = split("",$sseeqq);
for($iiii=1;$iiii<=scalar(@sdfgh);$iiii++){
print out "$sdfgh[$iiii]";
if($iiii%70==0) {print out "\n";}
}print out  "\n";
    # print out "$sseeqq\n";
            
    }   
  }
 }
}
