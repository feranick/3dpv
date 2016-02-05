#!/usr/bin/perl
# Converts .txt files of triangles (from the 3dpv code and misfitmodel format to Mathematica format
$txtfile = shift(@ARGV) ;
$output = "path-$txtfile" ;

open(TXT, "$txtfile") || die "Cannot open file $txtfile";
open(OUT, ">$output") || die "Cannot open file $output";

$j=0 ;

while (<TXT>) {
$j++;

#$point=$j + 90;
$point = $j + 4 ;
chomp;
$v1='';
$v2='';
$v3='';

($v1,$v2,$v3) = split(/\s+|,/,$_) ;
$scale = 2000000;

$p1 = $v1/$scale;
$p2 = $v2/$scale;
$p3 = $v3/$scale;

#print OUT "t" . "$j" . "=
#print OUT "t" ."$j = Graphics3D[Polygon[{{$v1,$v2,$v3},{$v4,$v5,$v6},{$v7,$v8,$v9}}]]\n" ;
print OUT "t" ."$point = Graphics3D[{PointSize[Large], Green, Point[{$p1, $p2, $p3}]}];\n" ;


}
#$point=$j+1 ;
#print OUT "t" ."$point = Graphics3D[{PointSize[Large], Red, Point[{5.0, 0.05, 5.0}]}]\n" ;

print OUT "Show[" ;
for ($k=1;$k <= $point-1; $k++) {
print OUT "t"."$k,";
}
print OUT "t$point]\n" ;
