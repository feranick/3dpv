#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# usage rotationz.pl filename degrees

$xyzfile = shift (@ARGV) ;
$final = "area.$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

#$side=10 ; ##should be used with cubic simulation box, or remove normalization of area 
$totarea = 0 ;
$j=0 ;
print FIN "aa={" ;
while (<XYZ>) {

$j++ ;
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';
          $v3='';
          $v4='';
          $v5='';
          $v6='';
          $v7='';
          $v8='';
          $v9='';
      
($v1,$v2,$v3,$v4,$v5,$v6,$v7,$v8,$v9) = split(/\s+|,/,$_) ;

$cross = sqrt(((-$v3)*$v5 + $v2*$v6 + $v3*$v8 - $v6*$v8 - $v2*$v9 + $v5*$v9)**2 + ($v3*$v4 - $v1*$v6 - $v3*$v7 + $v6*$v7 + $v1*$v9 - $v4*$v9)**2 +((-$v2)*$v4 + $v1*$v5 + $v2*$v7 - $v5*$v7 - $v1*$v8 + $v4*$v8)**2) ;
$area = ($cross/2);
$totarea = $totarea + $area ;
print FIN "$area," ; 
}
print FIN "}\n\n BinCounts[aa, {0, 100, 5}]\n" ;
print FIN "BarChart[%, BarSpacing -> 0]\n" ;
print FIN "\n The total area is $totarea m^2";




