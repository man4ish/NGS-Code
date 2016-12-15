#!/usr/bin/perl
open(fh,$ARGV[0]) or die "could not open the file\n";

while($input = <fh>)
{
   chomp($input);
     @a = split("\t", $input);
     if($a[8] !~/rs/) 
     {
        print $input,"\n";  
     }
}
