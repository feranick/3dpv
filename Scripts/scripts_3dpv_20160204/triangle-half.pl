#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# this algorithm doubles the number of triangles by halving them and redefining two from each of them
# if applied repeatedly, it can generate 2^n *t triangles, t=initial# triangles
# usage triangle-half.pl filename 

$xyzfile = shift (@ARGV) ;
$final = "d-"."$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

#TRANSLATION PART
$tx=0;
$ty=0;
$tz=0;

#SCALE FACTOR
$s=1;
$j=0;

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

          $p1x= ($v1 + $tx)*$s; 
          $p1y= ($v2 + $ty)*$s;
          $p1z= ($v3 + $tz)*$s;

          $p2x= ($v4 + $tx)*$s;  
          $p2y= ($v5 + $ty)*$s;
          $p2z= ($v6 + $tz)*$s;

          $p3x= ($v7 + $tx)*$s;  
          $p3y= ($v8 + $ty)*$s;
          $p3z= ($v9 + $tz)*$s;

#this subroutine doubles the triangle, and prints two from the original triangle
#define midpoint-vector
$hx = ($p2x - $p1x)/2;
$hy = ($p2y - $p1y)/2;
$hz = ($p2z - $p1z)/2; 
#define mid-point between P1 and P2
$mx = $p1x + $hx; 
$my = $p1y + $hy;
$mz = $p1z + $hz;

#define two triangles from the bisection
#triangle1 = s
$s1x = $p1x ; $s1y = $p1y ; $s1z = $p1z ; #first vertex
$s2x = $mx  ; $s2y = $my  ; $s2z = $mz  ; #second vertex
$s3x = $p3x ; $s3y = $p3y ; $s3z = $p3z ; #third vertex

#triangle2 = n
$n1x = $p2x ; $n1y = $p2y  ; $n1z = $p2z ; #first vertex
$n2x = $mx  ; $n2y = $my   ; $n2z = $mz  ; #second vertex
$n3x = $p3x ; $n3y = $p3y  ; $n3z = $p3z ; #third vertex

#Storing in new file: store as s1 s3 s2, n1 n3 n2 to perform "shuffling" at the same time"
# This will allow the next pass in the algorithm to bisect into smaller but similar triangles, instead of creating a pattern
# looking like a an array of rays.
printf FIN '%f', $s1x;
print FIN ",";
printf FIN '%f', $s1y;
print FIN ",";
printf FIN '%f', $s1z;
print FIN " ";

printf FIN '%f', $s3x;
print FIN ",";
printf FIN '%f', $s3y;
print FIN ",";
printf FIN '%f', $s3z;
print FIN " ";

printf FIN '%f', $s2x;
print FIN ",";
printf FIN '%f', $s2y;
print FIN ",";
printf FIN '%f', $s2z;
print FIN " ";
print FIN "\n" ;

printf FIN '%f', $n1x;
print FIN ",";
printf FIN '%f', $n1y;
print FIN ",";
printf FIN '%f', $n1z;
print FIN " ";

printf FIN '%f', $n3x;
print FIN ",";
printf FIN '%f', $n3y;
print FIN ",";
printf FIN '%f', $n3z;
print FIN " ";

printf FIN '%f', $n2x;
print FIN ",";
printf FIN '%f', $n2y;
print FIN ",";
printf FIN '%f', $n2z;
print FIN " ";
print FIN "\n" ;
}
