#!/bin/sh
. env.sh
export SYNASEARCH=/opt/synamatix/synasearch/synasearch-1.4.1-SNAPSHOT/bin/synasearch.sh
./SXGenMutatedSeq rna.fa transfile 
$SYNASEARCH -f transfile.mutatedseq.fna -mnl 15 -mrf 6.0 -fs 28 -fl 250 -fi 10 -fo 10 -mo 11 -Q -db hsr_gnom_ncbi-36.3_sb3.0.9 > snpseq.pw                                                         #choose the parameter  
./SXGenomeSnps snpseq.pw transfile.base 
