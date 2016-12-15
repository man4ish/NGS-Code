#! /usr/bin/env perl

###############################################################################
#
# fastqToFasta.pl
#
# Script to convert Solexa fastq file format, which includes both sequence
# tags and quality scores in "symbolic" format to a separate FASTA and QUAL
# file for easy loading through the existing Alpheus pipeline.
#
# (I'm holding off writing a real Java pipeline version until it's clearer what
# our data will usually look like)
#
# fastq file format:
#
#   @HWUSI-BETA8_1_1_1_947_886
#   AAGAAACCCTACATTTGTAGGGTTTGTCTTCAGTAT
#   +HWUSI-BETA8_1_1_1_947_886
#   YYYYXYYYXYYXXYYYYYXYWYYYXYYWYXUVRPPP
#
#  Where @ is header, followed by sequence, followed by quality header (+) 
#  followed by symbolic quality.
# 
#  From Irina at Solexa:
#  
#  "Quality scores are in symbolic format as a compact string of ASCII 
#   characters. Subtract 64 from the ASCII code to get the corresponding
#   quality values (in perl it's ord function), but bear in mind that under
#   the Solexa numbering scheme quality values can theoretically be as low
#   as -5, which has an ASCII encoding of 59=';'. for example 
#   ord( 'V' ) - 64 = 22 => the score of first base in all these
#   sequences is 22. etc...."
#
# Usage:  fastqToFasta.pl inputFile outputFile [run] [-n ] [ -v ]
#
# If supplied, the run argument will be added to FASTA header lines as
# "run=$run". 
#
# If supplied, the -n switch will cause all symbolic quality scores less 
# than zero to be normallized to zero. This is done purely to accommodate
# the Alpheus db which cannot currently handle negative values.
#
# If -v switch is supplied along with the -n switch, normalizing of negative
# values will be logged by printing to STDERR
#
# Will result in two files being written: outputFile.fna and outputFile.qual
# If run argument is supplied, run header info will be added to each header
# line as "run=$run".
# 
###############################################################################

use strict;
use Getopt::Long;

# Get the command line options and print a usage statement if not all required
# options are provided.
my $file;
my $outFile;
my $run;
my $normalize;
my $verbose;

my $result = GetOptions( "F=s" => \$file,
                         "O=s" => \$outFile,
                         "R=s" => \$run,
                         "n" => \$normalize,
                         "v" => \$verbose, );


if ( !$file || !$outFile ) {
    die( "\nUsage: fastqToFasta.pl -F inputFile -O outputFile [-R run] [-n] [-v]\n\n". 
         "    -F  The name of the input file in FASTQ format.\n".
         "    -O  The name of the output files.  Two output files are written\n".
         "                  named {outputFile}.fna and {outputFile}.qual.\n". 
         "    -R  If supplied, this argument represents the name of the run, and\n".
         "                  is added to the FASTA header as \"run=\$run\".\n".
         "    -n  If this switch is supplied, all symbolic quality scores less\n".
         "                  than zero are normalized to zero. This is done\n".
         "                  to accommodate the Alpheus DB, which cannot currently\n".
         "                  handle negative values.\n".
         "    -v  If this switch is supplied along with the -n switch, normalization\n".
         "                  of negative values is logged to STDERR.\n\n"
    );
}

my $fastaFile = "$outFile.fna";
my $qualFile  = "$outFile.qual";

open ( FILE, $file ) || die( "Error opening $file: $!\n" );
open ( FASTA, ">$fastaFile" ) || die( "Error opening $fastaFile: $!\n" );
open ( QUAL, ">$qualFile" ) || die( "Error opening $qualFile: $!\n" );

my @seqLines;

my $lineCounter = 0;
while( <FILE> ) {
    chomp;
    if ( $lineCounter == 4 ) {
        printSeqLines( \@seqLines, $run, $normalize, $verbose );
        
        $lineCounter = 0;
        undef( @seqLines );
    }

    push( @seqLines, $_ );
    $lineCounter++;
}

# print last entry from file
printSeqLines( \@seqLines, $run, $normalize, $verbose );

close( FILE ) || die( "Error closing $file: $!\n" );
close( FASTA ) || die( "Error closing $fastaFile: $!\n" );
close( QUAL ) || die( "Error closing $qualFile: $!\n" );

###############################################################################
#
# Print array of seq lines representing 4 lines of data from fastq file.
#
# printSeqLines( $arrayPtr, $run, $normalize, $verbose )
# 
# params: 
#        $arrayPtr  -> reference to seqLines fastq data array
#        $run       -> run label to append to FASTA headers (may be empty)
#        $normalize -> flags whether to normalize negative quality scores to 
#                      0 for use in alpheus
#        $verbose   -> flags whether to log normalizing of negative values
#
###############################################################################

sub printSeqLines( $$$ ) {
    my $arrPtr = shift;
    my @seqLines = @$arrPtr;
    my $run = shift;
    my $normalize = shift;
    my $verbose = shift;
            
    if ( $seqLines[ 0 ] =~ /^@(\S+)/ ) {
        my $header = $1;
        print FASTA ">$header";
        print QUAL ">$header";
        
        if ( $run ) {
            print FASTA " run=$run";
            print QUAL " run=$run";
        }
        
        print FASTA "\n";
        print QUAL "\n";

    } else {
        die( "Expected seq header line: $seqLines[ 0 ]\n" );
    }
    
    print FASTA "$seqLines[ 1 ]\n";
    
    my $first = 1;
    my @charArr = split( //, $seqLines[ 3 ] );
    for ( my $i = 0; $i < scalar( @charArr ); $i++ ) {
        my $convertedChar = ord( $charArr[ $i ] ) - 64;
        
        if ( $convertedChar < 0 && $normalize ) {
            # if verbose mode requested, print record of normalization
            if ( $verbose ) {
                print STDERR "Qual < 0: $convertedChar\n";
            }
            $convertedChar = 0;
        }
        
        if ( !$first ) {
            print QUAL " ";
        }
        print QUAL $convertedChar;
        $first = 0;
    }
    print QUAL "\n";
}
        