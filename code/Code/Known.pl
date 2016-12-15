#!/usr/bin/perl
open(fh,$ARGV[0]) or die "could not open the file\n";
 print "Chromosome\tGiNumber\tOffset\tNucleotide_Variant\tTotal_SNP_Reads\tTotal_Read_Density\tAvg_Read_Density\tAvg_QScore\tdbSNP\tGene_Name\tKeyword\tGene_Description\tExon\tZygosity\tLocal_CN\tKnown_CNVRegion\tTransversions(V)_Transitions(S)\tSynonymous(S)_Non-Synonymous(NS)\tProtein_Variant\tMissense\tNonsense\n";

while($input = <fh>)
{
   chomp($input);
     @a = split("\t", $input);
     if($a[8] =~/rs/) 
     {
        print $input,"\n";  
     }
}
