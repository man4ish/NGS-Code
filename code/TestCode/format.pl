#!/bin/sh
open(fh, $ARGV[0]) or die "could not open file\n";
$seq = "";
while($input = <fh>)
{
   chomp($input);
   if($input =~ />/)
   { 
     $seq.$input; 
   }
   print $seq,"\n";
}

