#!/usr/bin/perl
# This script rotates the coordinates of a 3dpv structure 
# usage rotationz.pl filename and input parameters by hand!

$xyzfile = shift (@ARGV) ;
$final = "roty"."$deg"."$xyzfile" ;

open(XYZ, "$xyzfile")  || die "Cannot open file $xyzfile";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

my $pi = 3.14159265358979;

sub deg_to_rad { ($_[0]/180) * $pi }
#ROTATION PART
#test OK
#print 'sin 30 degrees is ', sin(deg_to_rad(30)), "\n";
$theta = deg_to_rad($deg) ;
$r11= 1 ;
$r12= 0 ; 
$r13= 0 ;

$r21= 0 ; 
$r22= cos($theta) ; 
$r23= sin($theta) ;

$r31= 0;
$r32= -sin($theta) ; 
$r33= cos($theta) ;
#TRANSLATION PART
$ty=-8.700000;

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
$p1 = $v1*$r11 + $v2*$r12 + $v3*$r13 ;
$p2 = $v1*$r21 + $v2*$r22 + $v3*$r23 ;
$p3 = $v1*$r31 + $v2*$r32 + $v3*$r33 ;
$p4 = $v4*$r11 + $v5*$r12 + $v6*$r13 ;
$p5 = $v4*$r21 + $v5*$r22 + $v6*$r23 ;
$p6 = $v4*$r31 + $v5*$r32 + $v6*$r33 ;
$p7 = $v7*$r11 + $v8*$r12 + $v9*$r13 ;
$p8 = $v7*$r21 + $v8*$r22 + $v9*$r23 ;
$p9 = $v7*$r31 + $v8*$r32 + $v9*$r33 ;
## $result = sprintf("%08d", $number);
#$t1 = substr($p1, 0, 9);
#$t2 = substr($p2, 0, 9);
#$t3 = substr($p3, 0, 9);
#$t4 = substr($p4, 0, 9);
#$t5 = substr($p5, 0, 9);
#$t6 = substr($p6, 0, 9);
#$t7 = substr($p7, 0, 9);
#$t8 = substr($p8, 0, 9);
#$t9 = substr($p9, 0, 9);

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
#@transcoords[$j] = "$t1,$t2,$t3 $t4,$t5,$t6 $t7,$t8,$t9 \n"  ; 

}
#print FIN @transcoords; 
