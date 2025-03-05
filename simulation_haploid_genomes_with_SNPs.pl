#!/usr/bin/perl

use strict;
use warnings;
use File::Slurp;
use Getopt::Std;

# Function to display help message
sub usage {
    print <<"END_USAGE";
Usage: $0 -v <vcf_file> -r <ref_file> -s <sample_name> -o <output_prefix>

Options:
    -v <vcf_file>       Specify the input VCF file.
    -r <ref_file>       Specify the reference FASTA file.
    -s <sample_name>    Specify the sample name in the VCF file.
    -o <output_prefix>  Specify the output file prefix (e.g., 'hap' for SampleA_hap1.fa and SampleA_hap2.fa).
    -h                  Display this help message.
END_USAGE
    exit;
}

# Parse command-line options
my %options;
getopts('v:r:s:o:h', \%options);

# Display help if -h is provided or no arguments are given
usage() if $options{h} || !$options{v} || !$options{r} || !$options{s} || !$options{o};

# Input and output files
my $vcf_file = $options{v};       # VCF file
my $ref_file = $options{r};       # Reference FASTA file
my $sample_name = $options{s};    # Sample name in VCF file
my $output_prefix = $options{o};  # Output file prefix

# Read the reference FASTA file
my @ref_fasta = read_file($ref_file);
my %ref_seqs;
my $current_chrom;
for my $line (@ref_fasta) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
        $current_chrom = $1;
    } else {
        $ref_seqs{$current_chrom} .= $line;
    }
}

# Initialize haplotypes
my %hap1_seqs = %ref_seqs;  # hap1.fa starts with the reference sequences
my %hap2_seqs = %ref_seqs;  # hap2.fa starts with the reference sequences

# Open the VCF file
open(my $vcf_fh, '<', $vcf_file) or die "Cannot open VCF file $vcf_file: $!";

# Find the sample column index
my $sample_index = -1;
while (my $line = <$vcf_fh>) {
    chomp $line;
    if ($line =~ /^#CHROM/) {
        my @fields = split("\t", $line);  # 使用制表符分割
        for (my $i = 9; $i < @fields; $i++) {
            if ($fields[$i] eq $sample_name) {
                $sample_index = $i;
                last;
            }
        }
        last;
    }
}
die "Sample '$sample_name' not found in VCF file.\n" if $sample_index == -1;

# Process each variant in the VCF file
while (my $line = <$vcf_fh>) {
    chomp $line;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genotypes) = split("\t", $line);
    my $genotype = $genotypes[$sample_index - 9];
    my ($gt, @fields) = split(':', $genotype);

    # Replace bases in hap1.fa for 1/1 variants
    if ($gt eq '1/1') {
        substr($hap1_seqs{$chrom}, $pos - 1, length($ref)) = $alt;
    }

    # Replace bases in hap2.fa for 0/1 variants
    if ($gt eq '0/1') {
        substr($hap2_seqs{$chrom}, $pos - 1, length($ref)) = $alt;
    }
}
close($vcf_fh);

# Write hap1.fa
my $hap1_file = "${sample_name}_${output_prefix}1.fa";  # 输出文件名带上样本名
open(my $hap1_fh, '>', $hap1_file) or die "Cannot open file $hap1_file: $!";
for my $chrom (sort keys %hap1_seqs) {
    print $hap1_fh ">$chrom\n";
    while ($hap1_seqs{$chrom} =~ /(.{1,70})/g) {
        print $hap1_fh "$1\n";
    }
}
close($hap1_fh);

# Write hap2.fa
my $hap2_file = "${sample_name}_${output_prefix}2.fa";  # 输出文件名带上样本名
open(my $hap2_fh, '>', $hap2_file) or die "Cannot open file $hap2_file: $!";
for my $chrom (sort keys %hap2_seqs) {
    print $hap2_fh ">$chrom\n";
    while ($hap2_seqs{$chrom} =~ /(.{1,70})/g) {
        print $hap2_fh "$1\n";
    }
}
close($hap2_fh);

print "Haplotype files have been generated: $hap1_file and $hap2_file\n";
