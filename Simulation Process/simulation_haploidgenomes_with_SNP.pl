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
    -v <vcf_file>       Input VCF file containing SNPs
    -r <ref_file>       Reference genome FASTA file
    -s <sample_name>    Sample name in VCF (e.g., "AT01")
    -o <output_prefix>  Output prefix (e.g., "hap" generates SampleA_hap1.fa)
    -h                  Show this help message
END_USAGE
    exit;
}

# Parse command-line options
my %options;
getopts('v:r:s:o:h', \%options);

# Show help if requested or missing required args
usage() if $options{h} || !$options{v} || !$options{r} || !$options{s} || !$options{o};

# Input/Output files
my $vcf_file    = $options{v};      # Input VCF
my $ref_file    = $options{r};      # Reference FASTA
my $sample_name = $options{s};      # Sample name in VCF
my $out_prefix  = $options{o};      # Output prefix

# 1. Read reference genome
my @ref_fasta = read_file($ref_file);
my %ref_seqs;
my $current_chrom;
foreach my $line (@ref_fasta) {
    chomp $line;
    if ($line =~ /^>(\S+)/) {
        $current_chrom = $1;
    } else {
        $ref_seqs{$current_chrom} .= uc($line);  # Store as uppercase
    }
}

# 2. Initialize haplotypes
my %hap1 = %ref_seqs;  # Haplotype 1 (initialized as reference)
my %hap2 = %ref_seqs;  # Haplotype 2 (initialized as reference)

# 3. Process VCF file
open(my $vcf_fh, '<', $vcf_file) or die "Cannot open VCF file: $!";

# Find sample column index
my $sample_col = -1;
while (my $line = <$vcf_fh>) {
    next if $line !~ /^#CHROM/;
    my @headers = split("\t", $line);
    for (my $i = 9; $i < @headers; $i++) {
	$headers[$i] =~ s/[\"\'\s]//g;
        if ($headers[$i] eq $sample_name) {
            $sample_col = $i;
            last;
        }
    }
    last;
}
die "Sample '$sample_name' not found in VCF!" if $sample_col == -1;

# Process each variant
while (my $line = <$vcf_fh>) {
    next if $line =~ /^#/;  # Skip headers
    chomp $line;
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split("\t", $line);
    
    # Extract genotype (GT) for the target sample
    my $sample_gt = $samples[$sample_col - 9];
    my ($gt) = split(':', $sample_gt);  # Get GT (e.g., "0/1")
    
    # Skip if no genotype call (e.g., "./.")
    next if $gt !~ /^(\d+)[\/|](\d+)/;
    
    # Process SNPs (skip INDELs if length differs)
    next if length($ref) != length($alt);
    
    # Fix: 1/1 modifies BOTH haplotypes
    if ($gt eq '1/1') {
        substr($hap1{$chrom}, $pos - 1, length($ref)) = $alt;
        substr($hap2{$chrom}, $pos - 1, length($ref)) = $alt;
    } 
    # 0/1 modifies only haplotype 2
    elsif ($gt eq '0/1') {
        substr($hap2{$chrom}, $pos - 1, length($ref)) = $alt;
    }
}
close($vcf_fh);

# 4. Write output FASTA files
sub write_fasta {
    my ($file, $seqs) = @_;
    open(my $fh, '>', $file) or die "Cannot write to $file: $!";
    foreach my $chrom (sort keys %$seqs) {
        print $fh ">$chrom\n";
        my $seq = $$seqs{$chrom};
        while ($seq =~ /(.{1,80})/g) {  # 80 chars per line
            print $fh "$1\n";
        }
    }
    close($fh);
}

write_fasta("${sample_name}_${out_prefix}1.fa", \%hap1);
write_fasta("${sample_name}_${out_prefix}2.fa", \%hap2);

print "Generated:\n";
print "  - ${sample_name}_${out_prefix}1.fa (Haplotype 1)\n";
print "  - ${sample_name}_${out_prefix}2.fa (Haplotype 2)\n";
