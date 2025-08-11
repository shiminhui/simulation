#!usr/bin/perl
use strict;
use warnings;
#this perl script was used for filter the snp sites in which the total percent of missing data and low quality sample to all the population has reached a setting value.
#the low quality sample: if there is a sample has a PL value lower than 20, we identify a low quality genotype.
unless(@ARGV==5){
	die "Usage:\nperl $0 <snp_vcf_file> <output_file> <the_total_threshold_percent_of_bad_genotype_and_missing_data> <the_quality_score_for_filtering> <the sample number>\n"
}

my $vcf=shift;
my $out=shift;
my $badrate=shift;#the rate of missing data and low quality sample to the whole population, if the a snp site has a bad rate higher than this value, then this snp site will be filterred out.
my $qual=shift;
my $samplenum=shift;

open IN,"zcat $vcf|";
open OUT, "| gzip > $out";

while(my $line=<IN>){
	my $num0=0;
	my $num3=0;

	chomp($line);
	if($line=~/^#/){
		print OUT "$line\n";
		next;
	}#print header lines.
	my @temp=split /\t/,$line;

	foreach my $geno (9..$#temp){
		my $num1=0;
		my $num2=0;
		if($temp[$geno] =~ /^\.\/\./){#count the missing data number of each snp site.
			$num0++;
			next;
		}
		my @tmp=split /:/,$temp[$geno];
		my @pl=split /,/,$tmp[$#tmp];#split the three PL value of 0/0,0/1 and /1/1.
		foreach my $n (0..$#pl){
			if($pl[$n]==0){
				$num1++;
			}elsif($pl[$n]>0 and $pl[$n]<$qual){
				$num2++;
			}
		}
		if($num1==1 and $num2==0){#if there is a sample that has a PL with a 0 value and  two >20 values, then we say this is a qualified genotype.
                       $num3++;
                }
	}
	my $qualper=($samplenum-$num3)/$samplenum;#calculating the badrates (include missing data and bad genotype).
print "$temp[0]\t$temp[1]\t$qualper\n";	
	if($qualper<$badrate){
		print OUT "$line\n";
	}
}
close IN;
close OUT;
