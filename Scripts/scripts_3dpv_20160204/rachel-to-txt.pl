#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# usage rachel-to-txt.pl filename 

$xyzfile = shift (@ARGV) ;
$final = "conv"."$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

#TRANSLATION PART
$tx=0;
$ty=0;
$tz=0;

#SCALE FACTOR
$s=1;
$j=0 ;

while (<XYZ>) {
$j++ ;
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';
          $v3='';
	  $v4='';  # clear variables each time
          $v5='';
          $v6='';
          $v7='';  # clear variables each time
          $v8='';
          $v9='';
($v1,$v2,$v3,$v4,$v5,$v6,$v7,$v8,$v9) = split(/\s+|,/,$_) ;

          $p1= ($v1 + $tx)*$s; 
          $p2= ($v2 + $ty)*$s;
          $p3= ($v3 + $tz)*$s;

          $p4= ($v4 + $tx)*$s;  
          $p5= ($v5 + $ty)*$s;
          $p6= ($v6 + $tz)*$s;

          $p7= ($v7 + $tx)*$s;  
          $p8= ($v8 + $ty)*$s;
          $p9= ($v9 + $tz)*$s;

printf FIN '%f', $p1;
print FIN ",";
printf FIN '%f', $p2;
print FIN ",";
printf FIN '%f', $p3;
print FIN " ";

printf FIN '%f', $p4;
print FIN ",";
printf FIN '%f', $p5;
print FIN ",";
printf FIN '%f', $p6;
print FIN " ";

printf FIN '%f', $p7;
print FIN ",";
printf FIN '%f', $p8;
print FIN ",";
printf FIN '%f', $p9;
print FIN " ";
print FIN "\n" ;
}
