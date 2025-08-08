#!usr/bin/perl
use warnings;
use strict;

#unless(@ARGV==4){
#       die "perl $0 <i> <o> <low> <high>\n";
#}
#this perl script was used for filter the snp sites which have a dp value larger than  1100 and less than 6500.

my $in=shift;
my $out=shift;
my $low=shift;
my $high=shift;

open IN,"zcat $in|";
open OUT,"|/software/miniconda3/bin/bgzip > $out";

while(my $line=<IN>){
        chomp($line);
        my $number=0;
        print OUT "$line\n" if($line=~/^#/);
        if($line=~/;DP=(\d+);/){
                $number=$1;
                if($number > $low and $number < $high){
                        print OUT "$line\n";
                }
        }
}
close IN;
close OUT;
