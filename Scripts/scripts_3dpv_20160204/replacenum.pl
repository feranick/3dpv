#!/usr/bin/perl
# replaces e-06 with 
$txtfile = shift(@ARGV) ;
$output = "repl-$txtfile" ;
open(TXT, "$txtfile") || die "Cannot open file $txtfile";
open(OUT, ">$output") || die "Cannot open file $output";
$j=0 ;

while (<TXT>) {
chomp;
s/.e-06/0.000000/g;
print OUT "$_\n";
}
