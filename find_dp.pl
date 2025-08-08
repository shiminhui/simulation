#! /usr/bin/perl
use strict;
use warnings;

my $in=shift;
my $out=shift;
my $dpnum=shift;
open IN,"gzip -dc $in|";

open OUT,'>',$out;

my @aa;
my @dps;

while (my $line=<IN>) {
        chomp($line);
        my @line=split/\t/,$line;
        push @aa,$line[2];
}


@dps=sort {$a<=>$b} @aa;
print OUT "2.5%\t$dps[0.025*$dpnum-1]\n";
print OUT "97.5%\t$dps[0.975*$dpnum-1]\n";

close IN;
