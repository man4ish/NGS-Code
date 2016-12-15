#!/bin/sh
./SXSimulator_8Dec -f /home/jameswong/h_sapiens_v37_1/hs_grc_chr20.fa_masked -rand -o SIMrun20 -pe 205 5 -SNP 0.1 1 1 -s 30 -c 30 -dist
./SXSimulator_8Dec -f /home/jameswong/h_sapiens_v37_1/hs_grc_chr21.fa_masked -rand -o SIMrun21 -pe 205 5 -SNP 0.1 1 1 -s 30 -c 30 -dist
./SXSimulator_8Dec -f /home/jameswong/h_sapiens_v37_1/hs_grc_chr22.fa_masked -rand -o SIMrun22 -pe 205 5 -SNP 0.1 1 1 -s 30 -c 30 -dist
./SXSimulator_8Dec -f /home/jameswong/h_sapiens_v37_1/hs_grc_chrX.fa_masked -rand -o SIMrunX -pe 205 5 -SNP 0.1 1 1 -s 30 -c 30 -dist
./SXSimulator_8Dec -f /home/jameswong/h_sapiens_v37_1/hs_grc_chrY.fa_masked -rand -o SIMrunY -pe 205 5 -SNP 0.1 1 1 -s 30 -c 30 -dist
