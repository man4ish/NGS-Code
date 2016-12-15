#!/usr/bin/perl
open(fh, $ARGV[0]) or die "could not open the file\n";
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

%fields = {};
while($input= <fh>)
{
 @a = split(":",$input);
 $fields{$a[0]}= $a[1];
 #print "$a[0]\t$a[1]\n";
}

$dup = trim($fields{'Number of Duplicate Reads'});
$duppeecent = ""; 
if($dup=~/(\d+\.\d+)\%/)
{
$duppeecent = $1;
}

$tim = trim($fields{'Total Elaspse'});
$timetaken = "";
if($tim=~/(\d+\.\d+)/)
{
$timetaken = $1;
}

print  trim($fields{'File Name'}),"\t",trim($fields{'Number of Reads'}),"\t",trim($fields{'Average Length'}),"\t",trim($fields{'Total Base Pair'}),"\t",$timetaken,"\t",$duppeecent,"\t",trim($fields{'Minimum Quality Score'}),"\t",trim($fields{'Maximum Quality Score'}),"\t",trim($fields{'Minimum AVQS'}),"\t",trim($fields{'Maximum AVQS'}),"\t",trim($fields{'Bases QS >16   (%)'}),"\t",trim($fields{'Bases QS >20   (%)'}),"\t",trim($fields{'Reads with Average QS > 10 (%)'}),"\t",trim($fields{'Reads with Average QS > 16 (%)'}),"\t",trim($fields{'Reads with Average QS > 20 (%)'}),"\n";
