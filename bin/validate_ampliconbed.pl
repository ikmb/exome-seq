#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $infile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

chomp(my @lines = <$IN>);

my $lc = 0;

# Get all the information and perform value validation
foreach my $line (@lines) {

	$lc += 1;

	my @el = split("\t", $line);

        my $start = @el[1] ;
        my $stop = @el[2] ;
        my $len = $stop-$start ;

        $len < 40 or die "The length of the primer binding site in line $lc seems too large! Sure you are providing primer positions?\n"

	scalar(@el) == 6 or die "Incorrect number of elements detected in line $lc of your amplicon BED file (should be 6)!\n";

}

close($IN);

exit 0;
