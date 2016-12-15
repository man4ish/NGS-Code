#!/usr/bin/perl
sub manual{
   print "\nUsage:\n";
   print "~~~~~~~~~~~~~~~sRowFormatter~~~~~~~~~~~~~~~~~~\n";
   print "Perl RowFormatter.pl <input file> <sample type> >output file\n\n ";
   exit(0);
 }

sub match {
@qre1 = split('',$_[0]);
@qre2 = split('',$_[1]);
$scor=0;
for($w=0;$w<scalar(@qre1);$w++) {
  if($qre1[$w] eq $qre2[$w]) {
     $scor++;
       }  
    }
return $scor;
}

sub getquerypos
{
my $qen = $_[0];
$co =0;
@qendarr = split('',$qen);
for($u =0;$u<scalar(@qendarr);$u++)
{
  if($qendarr[$u] ne '-') {
      $co++;
   } 
 }
  return $co;
}

sub getlength
{
   $sum =0;
   $tar = $_[1];
   @target=split("",$tar);
      for($k =1;$k<=$_[0];$k++)   
      {
         if($target[$k] ne '-')
         {
            $sum += 1;
          }
       }
     return $sum;
 }

sub alignment {
  $q1=$_[0];
  $t1=$_[1];
  @q = split('',$q1);
  @t = split('',$t1);
  $start =$_[2];
  my  $pos=0; $delsum=0; $missum=0; $insum=0;
  @a="";
  if($t[0] eq '-')
  {
     $start = $start -1;
  }
  for($i=0;$i<scalar(@q);$i++) {
   if($q[$i] eq '-')       {  $delsum=0;  $deletionflag=0; $mismatchflag=0;  $insertionflag=1;
        $pos =$start + getlength($i,$t1);
         $y2 = "$pos"."."."-"."$t[$i]";
         push(@a,$y2);         
     }
   elsif($t[$i] eq '-') {
       if($deletionflag eq 1 && $deletionflag ne "")
        {$delsum +=1;
         }
        else {$delsum =0;}
        $deletionflag=1; $mismatchflag=0;$insertionflag=0;
        $r = getlength($i,$t1);
        $pos =$start + getlength($i,$t1);
        if($delsum!=0){ 
        $y3= "$pos"."."."$delsum"."$q[$i]"."-"; 
        push(@a,$y3);}
        else { $y4="$pos"."."."$q[$i]"."-";
             push(@a,$y4);
             }         
     }
   elsif($q[$i] ne $t[$i]) {$delsum=0; $deletionflag=0; $mismatchflag=1; $insertionflag=0;
        $pos =$start + getlength($i,$t1);
        $y6="$pos.$q[$i]$t[$i]";
        push(@a,$y6);        
     }
   else {$delsum=0; $deletionflag=0; $insertionflag=0; $mismatchflag=0;}
  }
 shift(@a);
 $count=scalar(@a);
$formatstr = "";
 for($id=0;$id<scalar(@a)-1;$id++){ 
    $formatstr .= $a[$id];
    $formatstr .=", ";
   }
  $formatstr .= $a[$count-1];
  if($formatstr ne "") {  print "\t$formatstr";}
}

if(scalar(@ARGV) eq 0) {
manual();
}
$count =0;
$filename=$ARGV[0];					
$outputfile=$ARGV[1];					
open(FILE,"$filename") or die "can not open file";
$samplename=$ARGV[1];					
chomp($samplename);
while($var=<FILE>){
	if($var=~/^sxoght:/){
		$str1=();$str2=();@ar=split(/\s+/,$var);
		if($ar[7]=~/\+/){
			$ar[7]=1;}
		else{
			$ar[7]=0;}
              if($ar[2] ne "hit"){
               print  ">$ar[2]\t$ar[8]\t$ar[9]\t$ar[1]\t";      
          }
	}	
	if($var=~/^Query:/){				
		$str1=$var;
		$str1=~s/^Query:\s+//g;
		$str1=~s/\d+\s+//g;
		$str1=~s/\s+//g;
	}
	if($var=~/^Sbjct:/){				
		$str2=$var;
		$str2=~s/^Sbjct:\s+//g;
		$str2=~s/\d+\s+//g;
		$str2=~s/\s+//g;
	}
	if(defined $str1 and defined $str2){
 $aa = getquerypos($str1);
 $matchscore = match($str1,$str2);
        if($ar[7] eq '1'){print "$ar[5]\t$aa\t$aa\t$ar[3]\t$ar[7]\t$ar[4]\t$matchscore";}
        else  {print "$aa\t$ar[5]\t$aa\t$ar[3]\t$ar[7]\t$ar[4]\t$matchscore"; }
           if (substr($str2,0,1) eq '-') {$ar = $ar[8]-1;}
		alignment($str1,$str2,$ar[8]);
	print  "\tT\tNM\n"; 
	($str1,$str2)=();$tmp=();
	}
}

