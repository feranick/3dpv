#!/usr/bin/perl
$framenumber = 0;
for ($n = 5; $n <= 9; $n++) {
for (0..4) {
$framenumber = $framenumber + 1;
   my $e = 2*$_;
  # print("Element $e was found at index $i\n");
#2010_6_15_8.61504.txt
$string = "perl color-scale.pl tower.txt 2011_6_30_$n.$e"."7801.txt $n.$e $framenumber";
print "\n$string";
`$string`;
}
}
for ($n = 10; $n <= 18; $n++) {
for (0..4) {
$framenumber = $framenumber + 1;
   my $e = 2*$_;
  # print("Element $e was found at index $i\n");
$string = "perl color-scale.pl tower.txt 2011_6_30_$n.$e"."78.txt $n.$e $framenumber";
print "\n$string";
`$string`;
}
}
