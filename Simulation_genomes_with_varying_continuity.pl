#!/usr/bin/perl

use strict;
use warnings;
use File::Slurp;
use Getopt::Std;

# Function to display help message
sub usage {
    print <<"END_USAGE";
Usage: $0 -i <input_file> -o <output_file>

Options:
    -i <input_file>    Specify the input FASTA file.
    -o <output_file>   Specify the output FASTA file.
    -h                Display this help message.
END_USAGE
    exit;
}

# Parse command-line options
my %options;
getopts('i:o:h', \%options);

# Display help if -h is provided or no arguments are given
usage() if $options{h} || !$options{i} || !$options{o};

# Input and output files
my $input_file = $options{i};
my $output_file = $options{o};

# Read the FASTA file into an array
my @fasta = read_file($input_file);

# Output the first line (header)
print $fasta[0];

# Merge the remaining lines into a single sequence string
my $merged_seq = join("", @fasta[1..$#fasta]);
$merged_seq =~ s/\s+//g;  # Remove all whitespace

# Find segments without 'N' and their positions
my @segments;
my $start = -1;
for (my $i = 0; $i < length($merged_seq); $i++) {
    my $char = substr($merged_seq, $i, 1);
    if ($char eq 'N') {
        if ($start >= 0) {
            push @segments, [$start, $i-1];
            $start = -1;
        }
    } else {
        $start = $i if $start < 0;
    }
}
if ($start >= 0) {
    push @segments, [$start, length($merged_seq)-1];
}

# Modify 300bp of normal sequence to 'N's at one-third of each segment
foreach my $seg (@segments) {
    my $len = $seg->[1] - $seg->[0] + 1;
    next if $len <= 900;  # Skip segments shorter than 900 bases
    my $pos = $seg->[0] + int($len / 3);
    $pos -= 1 if $len % 3 != 0;  # Round down if not divisible by 3
    # Replace 300bp of normal sequence with 'N's
    substr($merged_seq, $pos, 300) = "N" x 300;
}

# Open the output file and write the results
open(my $fh, '>', $output_file) or die "Cannot open file '$output_file': $!";
print $fh $fasta[0];     # Write the header
print $fh $merged_seq;   # Write the modified sequence
close($fh);

print "Results have been written to '$output_file'\n";
