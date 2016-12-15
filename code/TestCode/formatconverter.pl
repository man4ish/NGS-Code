#!/bin/perl
open(fh, $ARGV[0]) or die "could not open file\n";
$seq = "";
while($input = <fh>)
{
     if($input !~ />/)
     {
       chomp($input);
       $seq = $seq.$input;
     }
}
 print $seq,"\n";

