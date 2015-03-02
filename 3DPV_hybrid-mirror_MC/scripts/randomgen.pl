#!/usr/bin/perl
# Generates random set of triangles for 3dpv simulations ;
# Usage randomgen.pl namefile trinum xmax ymax zmax 
$txtfile = @ARGV[0] ;
$trinum = @ARGV[1] ;
$xmax= @ARGV[2] ;
$ymax= @ARGV[3] ;
$zmax= @ARGV[4] ;

for ($k=0, $k<5, $k++) {
  if (@ARGV[$k]==0) {

die "Usage: randomgen.pl namefile trinum xmax ymax zmax" ;

   }
 }
open(FIN, ">$txtfile") || die "Cannot open file $txtfile";

for ($j=1; $j < $trinum; $j++) {

$p1=rand($xmax);
 $p2=rand($ymax);
  $p3=rand($zmax);

$p4=rand($xmax);
 $p5=rand($ymax);
  $p6=rand($zmax);

$p7=rand($xmax);
 $p8=rand($ymax);
  $p9=rand($zmax);  

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
