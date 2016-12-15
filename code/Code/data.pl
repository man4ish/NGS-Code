
#!/usr/bin/perl
open(fh,$ARGV[0]) or "die could not open file\n";
while($input = <fh>)
{
  chomp($input);
  @a = split(/\t/, $input);
 if($a[76] ne $a[55])
{
print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$a[10]\t$a[11]\t$a[12]\t$a[13]\t$a[14]\t$a[15]\t$a[16]\t$a[17]\t$a[18]\t$a[19]\t$a[20]\t$a[25]\t$a[26]\t$a[27]\t$a[28]\t$a[34]\t$a[46]\t$a[47]\t$a[48]\t$a[49]\t$a[55]\t$a[67]\t$a[68]\t$a[69]\t$a[70]\t$a[76]\n";  
}
}
