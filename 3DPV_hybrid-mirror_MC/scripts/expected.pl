#!/usr/bin/perl
$powerfile = @ARGV[0] ;
$final = "expected.txt" ;
$j = 0;

open(POW, "$powerfile")  || die "Cannot open file $powerfile1";
#open(PPW, "$powerfile2")  || die "Cannot open file $powerfile2";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

$j = 0;
while (<POW>) { #sums over cells
$j++;
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';

($v1,$v2) = split(/\s+|,/,$_) ;
$e = $v2*0.96;
#$diff = @vector2[$j] - @vector1[$j] ;
print FIN "$v1  $e\n";
}
