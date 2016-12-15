#!/bin/sh
rm out
for i in `ls -a *.rank`
do
#echo $i
perl format.pl $i > $i.format
echo "./ParseSV " $i.format "CG1CG3 C SV >" $i.SVformat
done

