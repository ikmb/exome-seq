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
my $trio = undef;
my $help;

GetOptions(
    "help" => \$help,
    "trio=s" => \$trio,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

chomp(my @lines = <$IN>);

my $header = shift @lines;

my $ref_header = "patient;sample;status;library;readgroup;platform_unit;center;date;R1;R2";

my @elements = split(";" , $ref_header);

$header eq $ref_header or die "Samplesheet: header mal-formed!\n";

my $lc = 0;

my %valid_status= ( 0 => "true", 1 => "true" );

# Get all the information and perform value validation
foreach my $line (@lines) {

	$lc += 1;

	my @el = split(";", $line);

	scalar(@el) == scalar(@elements) or die "Incorrect number of elements detected in line $lc\n";

	my ($patient,$rgsm,$status,$rglb,$rgid,$rgpu,$center,$date,$r1,$r2) = @el ;

	if (!defined $valid_status{$status}) {
		die "Samplesheet: Found invalid status: $status \n";
	}
	
}

close($IN);

exit 0;
