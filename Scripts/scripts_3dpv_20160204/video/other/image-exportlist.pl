#!/usr/bin/perl
# Prepare list of images to be exported from Mathematica; change for particular case
$output = "savefiles.txt";
open(OUT, ">$output") || die "Cannot open file $output";
for($j=1; $j <= 70; $j++) {
print OUT "Export[\"~/work/3dpv/runs/Videos/ImagesVideo-tower/Image$j.tiff\",p$j];\n";
}
#print OUT "Export[\"~/work/3dpv/runs/Videos/ImagesVideo-tower/Image10.tiff\",p010];\n";
#$j=0;
#for ($j=11; $j< 71; $j++){
#$f = $j-10;
#print OUT "Export[\"~/work/3dpv/runs/Videos/ImagesVideo02/Image$j.tiff\",p$f];\n";
#}

