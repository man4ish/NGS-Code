#!/usr/bin/perl
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
for($u =0;$u<scalar(@qendarr);$u++) {
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
sub alignment
{
  $q1=$_[0];
  $t1=$_[1];
  @q = split('',$q1);
  @t = split('',$t1);
  $start =$_[2];
  my  $pos=0;
  $delsum=0;
  $missum=0;
  $insum=0;
  @a="";
  for($i=0;$i<scalar(@q);$i++)
  {
   if($q[$i] eq '-')       
     { 
       if($insertionflag eq 1 && $insertionflag ne "" )
        {
         $insum = $insum +1;            
         }
       else{$insum =0;}
        $delsum=0;
        $missum=0;
        $deletionflag=0;
        $mismatchflag=0;
        $insertionflag=1;
        $pos =$start + getlength($i,$t1);
        if($insum ne 0){ 
          $y1="$pos"."."."$insum"."-"."$t[$i]";
          push(@a,$y1);
           }
        else{ $y2 = "$pos"."."."-"."$t[$i]";
            push(@a,$y2);
             }         
     }
   elsif($t[$i] eq '-')
     {
       if($deletionflag eq 1 && $deletionflag ne "")
        {$delsum +=1;
         }
         else {$delsum =0;}
        $insum=0;
        $missum=0;
        $deletionflag=1;
        $mismatchflag=0;
        $insertionflag=0;
        $r = getlength($i,$t1);
        
        $pos =$start + getlength($i,$t1);
        if($delsum!=0){ 
        $y3= "$pos"."."."$delsum"."+"."$q[$i]"; 
           push(@a,$y3);}
        else { $y4="$pos"."."."+"."$q[$i]";
             push(@a,$y4);
             }         
     }
   elsif($q[$i] ne $t[$i])
     {
       if($mismatchflag ne 1 && $mismatchflag ne "")
        {$missum +=1;}
       else {$missum =0;}
        $insum=0;
        $delsum=0;
        $deletionflag=0;
        $mismatchflag=1;
        $insertionflag=0;
        $pos =$start + getlength($i,$t1);
       
        if($insum!=0){ $y5="$pos.$missum.$q[$i]$t[$i]";
            push(@a,$y5); }
        else { $y6="$pos.$q[$i]$t[$i]";
           push(@a,$y6);}         
     }
   else {$delsum=0;        
         $missum=0;
          $insum=0;
          $deletionflag=0;
          $insertionflag=0;
          $mismatchflag=0;}
  }
 shift(@a);
 $count=scalar(@a);
$formatstr = "";
 for($id=0;$id<scalar(@a)-1;$id++){ 
    $formatstr .= $a[$id];
    $formatstr .=", ";
   }
  $formatstr .= $a[$count-1];
  if($formatstr ne "") {
       print "\t$formatstr";
    }
}
$count =0;
sub manual{
   print "Manual RowFormatter.pl\n";
   print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
   print "Perl RowFormatter.pl <input file> <sample type> >output file\n\n ";
   exit(0);
 }
  open(fh,$ARGV[0]) or die "could not open the input file\n";

if(scalar(@ARGV) ne 2 )
{
  manual();
}
while($input=<fh>)
 {
  if($input=~/(\gi\|\d+\|ref.*\|)\s+((\-)?\d+)\s+(\d+\.\d+)\s+\d+\s+(\d+)\s+([-+])\s+\d+\s+\d+\s+(\d+)\s+(\d+)/) {
      push(@id,">$1");  
      push(@score,$2);
      push(@prb,$4);
      push(@match,$7);
      push(@orientaion,$6);

    }
  if($input=~/Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/) {
    push(@startarget,$1);
    push(@stoptarget,$3);
    push(@targetseq,$2);

   }
  if($input=~/Query:\s+(\d+)\s+(\S+)\s+(\d+)/) {
    push(@startquery,$1);
    push(@stopquery,$3);
    push(@queryseq,$2);
   }

  if($input=~/sxoght:\s+(\w\S+)\s/) {
     push(@illuminaid,$1);
    
   }
 }
$as=scalar(@id);
$bs=scalar(@targetseq);
$cs=scalar(@queryseq);
$ds=scalar(@illuminaid);

for($index =0;$index<scalar(@id);$index++) {
if($orientaion[$index] eq '+') {
print "$id[$index]\t$startarget[$index]\t$stoptarget[$index]\t$illuminaid[$index]\t$startquery[$index]\t";
print getquerypos($queryseq[$index]);
print "\t",getquerypos($queryseq[$index]);
print "\t$score[$index]\t";
print "1\t";
}
else {
print "$id[$index]\t$startarget[$index]\t$stoptarget[$index]\t$illuminaid[$index]\t";
print getquerypos($queryseq[$index]);
print "\t$startquery[$index]\t";
print getquerypos($queryseq[$index]);
print "\t$score[$index]\t";
print "0\t";
}
print "$prb[$index]\t";
print match($queryseq[$index],$targetseq[$index]);
alignment($queryseq[$index],$targetseq[$index],$startarget[$index]); 
print"\tO\t$ARGV[1]\n";
}

