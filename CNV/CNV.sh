#!/bin/sh
rm out
for i in `ls -a *.Annot`
do
echo $i
perl parse.pl $i CG4 N CNV > $i.format
done

