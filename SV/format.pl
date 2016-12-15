#!/usr/bin/perl
open (fh , $ARGV[0]) or die "could not open file\n";
while($input = <fh>)
{
  @a = split("\t",$input);
  $d = @a;
# if($a[0] ne "#ChrId1")
# { 
  for ($i = 0; $i < $d-1; $i++ )
  {
    if($a[$i] eq "")
    {
      print "-\t";
    }  
    else 
    {
      print "$a[$i]\t";
    }
  }
    chomp($a[$d-1]);
    if($a[$d-1] eq "")
    {
     print "-\n";
    }
     else
    {
      print "$a[$d-1]\n";
    }
# }
}
