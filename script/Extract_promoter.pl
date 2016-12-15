#!/usr/bin/perl -w
use strict;
use LWP::Simple;
############################################################################################################
#  perl Extract_promoter.pl -chr 1 -l 43143512351 -r 199 > promodb                                         #
#  Output:                                                                                                 # 
#  >ADAM30:chr1, prom_start= 119583635, prom_end= 119584635, Prom_orient= -                                #
#  GTGTATAAATGTTGAAACCATCCCTGATTTGCCAGAGCATACGACTATAA                                                      #
############################################################################################################
my $chr_num = $ARGV[1];
my $chr_length = $ARGV[3];
my $range = $ARGV[5];
my $argc = scalar(@ARGV);
sub manual
{   print "\nExtract_promoter.pl\n";
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print "perl Extract_promoter.pl -chr 1 -l 43143512351 -r 199 > promodb\n\n"; 
}

sub Getseq
{
    my $hdr = "http://rulai.cshl.edu/cgi-bin/TRED/tred.cgi?process=promInfo&pid=";
    my $gene_name = "";  
    my $catalog = get($hdr.$_[0]);
    die "Couldn't get it?!" unless $catalog;
    while ($catalog =~ /geneInfo&gid=.*>(.*)<.*a><.*td>/gi)
    {
           $gene_name = $1;
    }
    print ">$gene_name:chr$chr_num, prom_start= $_[1], prom_end= $_[2], chr_orient= $_[3]\n";
    while ($catalog =~ /"black">(.*)<.*font>&/gi)
    {
           print $1,"\n";
    }
}

sub Getindex
{
    my @genename =();
    my @orientation=();
    my @start=();
    my @promoterid=(); 
    my $info = "http://rulai.cshl.edu/cgi-bin/TRED/tred.cgi?process=searchProm&sel_organism=human&sel_type=genomic_coordinate&tx_search_terms=&chrom=";
    my $from = "&from=1&to=";
    my $rage = "&start=";
    my $ur = "&anotherPage=GO";  
    my $catalog = get($info.$chr_num.$from.$_[0].$chr_length.$range.$ur);
    die "Couldn't get it?!" unless $catalog;
    while ($catalog =~ /("row0">|"row1">)(A.*)<.*td>/gi)               #getting genename
    {
           push(@genename,$2);
    }
    while ($catalog =~ /("row0">|"row1">)(\+|\-)<.*td>/gi)             #getting orientation
    {
           push(@orientation,$2);
    }
    while ($catalog =~ /("row0">|"row1">)(\d+)<.*td>/gi)               #getting start position
    {
           push(@start,$2);
    }
    while ($catalog =~ /promInfo&pid=(.*)">/gi)                        #getting promoterid
    {         
          push(@promoterid,$1);
    }
    for(my $i =0; $i< scalar(@promoterid);$i++)
    {
        my $prm_end="";
        my $prm_start = "";
        if($orientation[$i] eq '-')
        {
            $prm_start = $start[$i]+700;
            $prm_end  =  $start[$i]-299;
        }
        elsif($orientation[$i] eq '+')
        {
            $prm_start = $start[$i]-700;
            $prm_end  =  $start[$i]+299;
        }
        Getseq($promoterid[$i],$prm_start,$prm_end,$orientation[$i]);
        print "\n";         
    }
}
 
if ($argc != 6)
{
   manual();
   exit();
}
for (my $k =0;$k<=$range; $k++)
{
     Getindex(20*$k);
}
