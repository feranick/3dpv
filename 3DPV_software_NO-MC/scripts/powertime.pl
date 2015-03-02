#!/usr/bin/perl
$output = "power-vs-time.txt";
open(OUT, ">$output") || die "Cannot open file $output";

#define subroutine to sum power from each file and return it into variable $sum
#pass time and name of power file

sub sumpower {
$sum = 0;
#$powfile =  ;
open(POW, "$powfile") || die "Cannot open file $txtfile";

while (<POW>) {
($p1,$p2) = split(/\s+|,/,$_) ;
$sum = $p2 + $sum;
}
$sum = $sum + 0; #ensure output right for now
}

# MAIN # 
# Need to change the part containing the filename
# as many for cycles as needed to cover all files
for ($hour = 6; $hour <= 9; $hour++) {
for ($min = 0; $min <= 4; $min++) {
   my $d = 2*$min  ;
$powfile = "2010_9_19_$hour.$d"."7958.txt";
$time = "$hour.$d" ;
$c = &sumpower;
print OUT "$time $c\n"; # is this enough to call the subroutine
}
}

for ($hour = 10; $hour <= 16; $hour++) {
for ($min = 0; $min <= 4; $min++) {
   my $d = 2*$min ;
$powfile = "2010_9_19_$hour.$d"."796.txt";
$time = "$hour.$d" ;
$c = &sumpower;
print OUT "$time $c\n"; # is this enough to call the subroutine
}
}
