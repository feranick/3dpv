#!/usr/bin/perl
$powerfile1 = @ARGV[0] ;
$powerfile2 = @ARGV[1] ;
$final = "difference.txt" ;
$j = 0;

open(POW, "$powerfile1")  || die "Cannot open file $powerfile1";
open(PPW, "$powerfile2")  || die "Cannot open file $powerfile2";
open(FIN, ">$final") || die "Cannot open file $final"; #APPENDS TEXT AT THE END OF IT

#first argument is the file with smaller power
while (<POW>) { #sums over cells
$j++;
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';
($v1,$v2) = split(/\s+|,/,$_) ;
@vector1[$j] = $v2;

}
$j = 0;
while (<PPW>) { #sums over cells
$j++;
chomp;    #remove trailing whitespaces

          s/^\s+//;  #remove leading whitespaces
          $v1='';  # clear variables each time
          $v2='';
($v1,$v2) = split(/\s+|,/,$_) ;
@vector2[$j] = $v2;
$diff = @vector2[$j] - @vector1[$j] ;
print FIN "$v1  $diff\n";
}
