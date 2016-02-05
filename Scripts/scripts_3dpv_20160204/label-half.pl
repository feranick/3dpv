#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# usage rachel-to-txt.pl filename 

$xyzfile = shift (@ARGV) ;
$final = "d-"."$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

while (<XYZ>) {
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';
($v1,$v2) = split(/\s+|,/,$_) ;
#this line doubles the labels like it's done for the files
print FIN "$v1 $v2 \n$v1 $v2 \n";
}
