#!/usr/bin/perl
$filename = "2011_4_21_";
$structfile = "Fl_spr_n1p4.txt";
$folderimages = "~/Desktop/test";
$framenumber = 0;
$n_var_in = 5;
$n_nvar_fin = 18;
$lab1 = "0839";
$lab2 = "084";
$numframes = 0;

for ($n = $n_var_in; $n <= 9; $n++) {
    for (0..4) {
        $framenumber = $framenumber + 1;
        $numframes = $numframes + 1;
        my $e = 2*$_;
        # print("Element $e was found at index $i\n");
        #2010_6_15_8.61504.txt
        #$string = "perl color-scale.pl cube35mm.txt 2011_11_24_$n.$e"."6162.txt $n.$e $framenumber";
        $string = "perl color-scale.pl $structfile $filename$n.$e"."$lab1.txt $n.$e $framenumber";
    
        print "\n$string";
        `$string`;
        }
    }
for ($n = 10; $n <= $n_nvar_fin; $n++) {
    for (0..4) {
        $framenumber = $framenumber + 1;
        $numframes = $numframes + 1;
        my $e = 2*$_;
        # print("Element $e was found at index $i\n");
        #$string = "perl color-scale.pl cube35mm.txt 2011_11_24_$n.$e"."616.txt $n.$e $framenumber";
        $string = "perl color-scale.pl $structfile $filename$n.$e"."$lab2.txt $n.$e $framenumber";
        print "\n$string";
        `$string`;
        }
    }




# Prepare list of images to be exported from Mathematica; change for particular case
$output = "dumpall.txt";
open(OUT, ">>$output") || die "Cannot open file $output";

print OUT "list={";
#for ($m = 1; $m <= 10; $m++) {
#   print OUT "p00$m,";
#   print OUT "p00$m,";
#   if ($m == 10) {
#       print OUT "p0$m,";
#   }
#}

for ($n = 1; $n <= $numframes; $n++) {
    print OUT "p$n,";
    print OUT "p$n,";
}
print OUT "}\n\n";


#$output = "dumpall.txt";
#open(OUT, ">>$output") || die "Cannot open file $output";
#print OUT "video = ListAnimate[list, 1, FrameMargins -> Large]\n";
#print OUT "Export[\"3DPVideo.swf\", list]";


$output = "dumpall.txt";
open(OUT, ">>$output") || die "Cannot open file $output";
for($j=1; $j <= 70; $j++) {
    print OUT "Export[\"$folderimages/Image$j.tiff\",p$j];\n";
}

