#!/usr/bin/perl 
open(fh,$ARGV[0]) or die "could not open the file\n";
%hash=();
while($input = <fh>) {
chomp($input);
$hash{$input}= "";
}

while ( ($k,$v) = each %hash ) {
    print "$k\n";
}
