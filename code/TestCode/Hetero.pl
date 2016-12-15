#!/usr/bin/perl
open(fh,$ARGV[0]) or die "could not open the file\n";
print "Chromosome\tGiNumber\tOffset\tNucleotide_Variant\tSNP_Reads\tTotal_Read_Density\tAvg_QScore\tdbSNP\tGene_Name\tGene_Description\tKeyword\tmiRNA\tPromoter\tExon\tZygosity\tLocal Copy Number\tKnown_CNVRegion\tTransversions(Tv)_Transitions(Ts\tSynonymous(S)_Non-Synonymous(NS)\tProtein_Variant\tMissense\tNonsense\n"; 
while($input = <fh>)
{
   chomp($input);
     @a = split("\t", $input);
     if($a[14] eq "het") 
     {
        print $input,"\n";  
     }
}
