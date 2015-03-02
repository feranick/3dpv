#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# usage rotationz.pl filename degrees
use Math::Trig;
#use Math::Trig 'deg2rad';
$xyzfile = shift (@ARGV) ;
$final = "arrows.$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(OUT, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

#set center of sun coordinate
#$xp= 4760444.172256 ;
#$yp= -48532275.304197  ;
#$zp= 140794021.268946 ;
#$l = sqrt($xp**2 + $yp**2 + $zp**2);

#And make it a unit vector
#$xpnorm = $xp / $l ;
#$ypnorm = $yp / $l ;
#$zpnorm = $zp / $l ;

#$totarea = 0 ;

$j=0 ;

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


#list components normal unit vector
$n1 = ((-$v3)*$v5 + $v2*$v6 + $v3*$v8 - $v6*$v8 - $v2*$v9 + $v5*$v9) / $cross ;
$n2 = ($v3*$v4 - $v1*$v6 - $v3*$v7 + $v6*$v7 + $v1*$v9 - $v4*$v9) / $cross ;
$n3 = ((-$v2)*$v4 + $v1*$v5 + $v2*$v7 - $v5*$v7 - $v1*$v8 + $v4*$v8) / $cross ;

$cm1 = ($v1 + $v4 + $v7)/3 ;
$cm2 = ($v2 + $v5 + $v8)/3 ;
$cm3 = ($v3 + $v6 + $v9)/3 ;

$fin1 = 3*$n1;
$fin2 = 3*$n2;
$fin3 = 3*$n3;

##Good line, to be restored after test
#print OUT "t" ."$j = Graphics3D[{Blue,Thick,Arrowheads[Medium],Arrow[{{$cm1,$cm2,$cm3},{$fin1,$fin2,$fin3}}], Boxed->False}]\n" ;
print OUT "t" ."$j = Graphics3D[{Blue,Thick,Arrowheads[Medium],Arrow[{{0,0,0},{$n1,$n2,0}}]}]\n" ;
}

$point = $j+1 ;
$cuboid = $j+2;
print OUT "t" ."$point = Graphics3D[{Opacity[0.2], Sphere[{0, 0, 0}], Boxed->False}]\n";
#print OUT "t" ."$cuboid = Graphics3D[{Opacity[0], Cuboid[{0, 0, 0}, {10, 10, 10}]}, Boxed->False]\n";
print OUT "t" ."$cuboid = Graphics3D[{PointSize[Large], Red, Point[{0.0, 0.0, 1.0}]}]\n" ;
##################
#Write Show line
print OUT "Show[" ;

for ($k=1;$k <= $j-1; $k++) {
print OUT "t"."$k,";
}
#print OUT "t$j,t$point,t$cuboid,Boxed->False]\n" ;
print OUT "t$j,t$point,t$cuboid,Boxed->False,ViewPoint->{0,0,1}]\n" ;
