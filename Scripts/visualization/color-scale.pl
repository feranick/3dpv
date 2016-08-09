#!/usr/bin/perl
# Usage: perl ./color-scale.pl structurefile.txt powerfile.txt time framenumber pov_x pov_y pov_z
#use 4-digit time format 
use List::Util qw(min max);
#$min = min @list;
 
$txtfile = @ARGV[0] ;
$powfile = @ARGV[1] ;
$time = @ARGV[2] ;
$framenumber = @ARGV[3] ;
$pov_x = @ARGV[4] ;
$pov_y = @ARGV[5] ;
$pov_z = @ARGV[6] ;

#$output = "$time-$txtfile" ;
$output = "dumpall.txt";
open(TXT, "$txtfile") || die "Cannot open file $txtfile";
open(POW, "$powfile") || die "Cannot open file $txtfile";
open(OUT, ">>$output") || die "Cannot open file $output";

##DEFINE MAX POWER SO TO RESCALE ALL COLORS ACCORDINGLY
#Color usually normalized at noon highest power 
$maxpower = 140 ;

#Cycle through power file and assign fractional colors to cells
$j=0 ; $s=0;
while (<POW>) {
$s++ ;
($p1,$p2) = split(/\s+|,/,$_) ;
if ($p2 < 0.001) { $p2 = 0; } #avoids expon form
$r=$p2/$maxpower; if ($r < 0.001) {$r=0; } #avoids expon form
@power[$s] = min(1,$r) ; #guarantees f < 1
}

while (<TXT>) {
$j++;

chomp;
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

#OKOK#print OUT "t" ."$j = Graphics3D[{Opacity[@power[$j]],Red,Polygon[{{$v1,$v2,$v3},{$v4,$v5,$v6},{$v7,$v8,$v9}}],Boxed->False}];\n" ;
print OUT "t" ."$j = Graphics3D[{Opacity[@power[$j]],Red,{EdgeForm[Thick],Polygon[{{$v1,$v2,$v3},{$v4,$v5,$v6},{$v7,$v8,$v9}}]}},Boxed->False,PlotLabel->$powfile];\n" ;
#print OUT "t" ."$j = Graphics3D[{Darker[Red, @power[$j]],Polygon[{{$v1,$v2,$v3},{$v4,$v5,$v6},{$v7,$v8,$v9}}],Boxed->False}]\n" ;

#print OUT "t" ."$point = Graphics3D[{Opacity[0.2], Sphere[{0, 0, 0}], Boxed->False}]\n";
}
$point=$j + 1 ;
$cuboid=$point + 2;
#print OUT "t" ."$point = Graphics3D[{PointSize[Large], Red, Point[{5.0, 0.05, 5.0}]}];\n" ;
#print OUT "t" ."$cuboid = Graphics3D[{Opacity[0], Cuboid[{0, 0, 0}, {10, 10, 10}]}, Boxed->False];\n";
print OUT "p$framenumber=Show[" ;
for ($k=1;$k <= $j-1; $k++) {
print OUT "t"."$k,";
}

#print OUT "t$j,t$point,t$cuboid,ViewPoint -> {1, -1, 1}]\n" ;
print OUT "t$j,ViewPoint -> {$pov_x, $pov_y, $pov_z}];\n\n" ;
