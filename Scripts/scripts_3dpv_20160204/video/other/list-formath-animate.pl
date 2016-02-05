#!/usr/bin/perl
$output = "list.txt";
open(OUT, ">$output") || die "Cannot open file $output";

print OUT "list={";
for ($m = 1; $m <= 10; $m++) {
print OUT "p00$m,";
print OUT "p00$m,";
if ($m == 10) {
print OUT "p0$m,";
print OUT "p0$m,";
}
}
for ($n = 1; $n <= 60; $n++) {
print OUT "p$n,";
print OUT "p$n,";
}
print OUT "}";
