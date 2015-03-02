#!/usr/bin/perl
# Converts .txt files of triangles (from the 3dpv code and misfitmodel format to Mathematica format
$txtfile = shift(@ARGV) ;
$output = "math-$txtfile" ;

open(TXT, "$txtfile") || die "Cannot open file $txtfile";
open(OUT, ">$output") || die "Cannot open file $output";

$j=0 ;

while (<TXT>) {
$j++;

chomp;
s/.e-06/0.000000/g;
$v1='';
$v2='';
$v3='';
$v4='';
$v5='';
$v6='';
$v7='';
$v8='';
$v9='';

($v1,$v2,$v3,$v4,$v5,$v6,$v7,$v8,$v9) = split(/\s+|,/,$_) ;
#print OUT "t" . "$j" . "=
print OUT "t" ."$j = Graphics3D[{EdgeForm[{Thick,Black}],Polygon[{{$v1,$v2,$v3},{$v4,$v5,$v6},{$v7,$v8,$v9}}]}];\n" ;
}
$point=$j+1 ;
print OUT "t" ."$point = Graphics3D[{PointSize[Large], Red, Point[{5.0, 0.05, 5.0}]}];\n" ;
print OUT "Show[" ;
for ($k=1;$k <= $j-1; $k++) {
print OUT "t"."$k,";
}
print OUT "t$j,t$point]\n" ;
