#!/usr/bin/perl
open(fh,$ARGV[0]) or die "could not open file\n";
while($input = <fh>)
{
  chomp($input);
  @a = split("\t",$input);
  print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$a[10]\t$a[11]\t$a[12]\t$a[15]\t$a[17]\n";
}

