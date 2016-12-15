#!/usr/bin/perl
open(fh,$ARGV[0]) or die "could not open the file\n";
while($input = <fh>)
{
   chomp($input);
   if($input !~ /#/ && (length($input) != 0) && $input !~/Chromosome/)
   {
     @a = split("\t", $input);
     if($a[12] eq "-")
     {
        print $input,"\n";
     }
   }
}
