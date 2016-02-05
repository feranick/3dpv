#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# usage rotationz.pl filename degrees
use Math::Trig;
#use Math::Trig 'deg2rad';
$xyzfile = shift (@ARGV) ;
$final = "hist.$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

#set center of sun coordinate
$xp= -92624559.406358 ;
$yp= -109365574.653009   ;
$zp= 40753675.608406  ;
$l = sqrt($xp**2 + $yp**2 + $zp**2);

#And make it a unit vector
$xpnorm = $xp / $l ;  
$ypnorm = $yp / $l ;
$zpnorm = $zp / $l ;

#$totarea = 0 ;

$j=0 ;
print FIN "BinCounts[{" ;
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

#test: print FIN "$n1 $n2 $n3\n";
#calculate angle between P vector and n vector:
$angle = acos($n1*$xpnorm + $n2*$ypnorm + $n3*$zpnorm) ;
$radangle = rad2deg($angle) ;

if ($radangle > 180) {
$radangle = $radangle - 180 ;
}

if ($radangle > 90) {
$radangle = 180 - $radangle;
}

#print FIN "{$radangle,$area},"
print FIN "$radangle,"
}
print FIN "},{0,90,30}]\n" ;
print FIN "BarChart[%, BarSpacing -> 0]\n\n";





