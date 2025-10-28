#!/usr/bin/perl

use strict;
use warnings;

# Configuration parameters
my $min_anchor_size = 4000000;  # Minimum anchor size for merging (4Mb)
my $max_gap_size    = 500000;   # Maximum allowed gap size between regions (500kb)

# Input validation
my $input_file = $ARGV[0];
die "Usage: perl $0 <input.bed>\n" unless defined $input_file;

# Read and parse input BED file
my @regions;
open my $fh, '<', $input_file or die "Cannot open file '$input_file': $!";
while (my $line = <$fh>) {
    chomp $line;
    # Skip empty lines and header lines
    next if $line =~ /^\s*$/ or $line =~ /^track/ or $line =~ /^#/; 
    
    my @fields = split /\t/, $line;
    die "Invalid format: line has fewer than 3 columns: $line" if @fields < 3;
    push @regions, \@fields;
}
close $fh;

# Exit if no valid regions found
exit 0 if !@regions;

# Sort regions by chromosome and start position
@regions = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @regions;

# Merge adjacent regions algorithm
my $merged_in_pass = 1;  # Flag to track if merging occurred in current pass
while ($merged_in_pass) {
    $merged_in_pass = 0;  # Reset flag for new pass
    
    my $i = 0;
    while ($i < $#regions) {  # Process all region pairs
        my $current = $regions[$i];
        my $next    = $regions[$i + 1];

        # Check if regions are on same chromosome and gap is below threshold
        if (
            $current->[0] eq $next->[0] and
            ($next->[1] - $current->[2]) < $max_gap_size
        ) {
            my $current_len = $current->[2] - $current->[1];
            my $next_len    = $next->[2] - $next->[1];

            # Merge if either region meets minimum anchor size
            if (($current_len > $min_anchor_size) or ($next_len > $min_anchor_size)) {
                $current->[2] = $next->[2];          # Extend current region
                splice(@regions, $i + 1, 1);         # Remove merged region
                $merged_in_pass = 1;                 # Flag that merging occurred
                next;                                # Skip increment to check new neighbor
            }
        }
        $i++;
    }
}

# Output merged regions
foreach my $region (@regions) {
    print join("\t", @$region), "\n";
}
