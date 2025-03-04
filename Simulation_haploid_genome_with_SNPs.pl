#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

# Function to display help message
sub usage {
    print <<"END_USAGE";
Usage: $0 -s <snp_file> -f <fasta_file> -o <output_file> -e <error_file>

Options:
    -s <snp_file>     Specify the SNP file (format: header\\tposition\\toriginal_base\\tnew_base).
    -f <fasta_file>   Specify the input FASTA file.
    -o <output_file>  Specify the output FASTA file after base replacement.
    -e <error_file>   Specify the error log file for mismatched bases.
    -h                Display this help message.
END_USAGE
    exit;
}

# Parse command-line options
my %options;
getopts('s:f:o:e:h', \%options);

# Display help if -h is provided or no arguments are given
usage() if $options{h} || !$options{s} || !$options{f} || !$options{o} || !$options{e};

# Input and output files
my $snp_file = $options{s};      # SNP file
my $fasta_file = $options{f};    # Input FASTA file
my $output_file = $options{o};   # Output FASTA file
my $error_file = $options{e};    # Error log file

# Read SNP file and store SNP data
open(my $snp_fh, "<", $snp_file) or die "Cannot open SNP file $snp_file: $!";

my %snps;
while (my $line = <$snp_fh>) {
    chomp $line;
    my ($header, $pos, $original_base, $new_base) = split("\t", $line);
    $snps{$pos} = {
        original_base => $original_base,
        new_base      => $new_base,
        header        => $header
    };
}
close $snp_fh;

# Read input FASTA file and store sequence data
open(my $fasta_fh, "<", $fasta_file) or die "Cannot open FASTA file $fasta_file: $!";

my $header = <$fasta_fh>;
chomp $header;
my $sequence = '';
while (my $line = <$fasta_fh>) {
    chomp $line;
    $sequence .= $line;
}
close $fasta_fh;

# Replace bases and log errors
open(my $error_fh, ">", $error_file) or die "Cannot open error file $error_file: $!";
foreach my $pos (sort {$a <=> $b} keys %snps) {
    my $original_base = $snps{$pos}->{original_base};
    my $new_base = $snps{$pos}->{new_base};
    my $header = $snps{$pos}->{header};

    if (substr($sequence, $pos - 1, 1) eq $original_base) {
        substr($sequence, $pos - 1, 1) = $new_base;
    } else {
        print $error_fh "$header\n";
    }
}
close $error_fh;

# Write the modified sequence to the output FASTA file
open(my $output_fh, ">", $output_file) or die "Cannot open output file $output_file: $!";
print $output_fh "$header\n";
while ($sequence =~ /(.{1,70})/g) {
    print $output_fh "$1\n";
}
close $output_fh;

print "Base replacement completed! Modified sequence saved to $output_file. Errors logged to $error_file.\n";
