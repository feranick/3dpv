#!/usr/bin/perl
$powerfile = shift (@ARGV) ;
$final = "sumpower_$powerfile" ;
$sum = 0;
open(POW, "$powerfile")  || die "Cannot open file $powerfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT
while (<POW>) { #sums over cells
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';
($v1,$v2) = split(/\s+|,/,$_) ;
$v2 = ($v2*10)/6; #rescaling to 10% efficiency ###### NOTA BENE!! ######
$sum = $sum + $v2;
}
print FIN "\n\n Total power at this time is $sum Watts\n\n";
print FIN "\n\n Total power at this time is $sum Watts\n\n";
