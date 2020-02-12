#!/bin/env perl

use Cwd;
use Getopt::Long;

my $usage = qq{
  Getting help:
    [--help] 
    [--region]
	A BED file with the target region to check

  Ouput:
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $region = undef;
my $help;

GetOptions(
    "help" => \$help,
    "region=s" => \$region,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my $dir = getcwd;

my $header = qq(
id: 'Coverage of the SRY gene'\n
section_name: 'ikmb/exome-seq sex check'\n
plot_type: 'html'\n
description: ' shows the coverage of thr SRY gene to determine sex of the sample.'\n
data: |\n
  <dl class="dl-horizontal">
);

printf $header . "\n";

opendir(DIR, $dir) or die $!;

while (my $file = readdir(DIR) ) {

	next unless ($file =~ m/\.bam$/);
	
	my $sample = (split /\./, $file )[0] ;

	my $data = `samtools bedcov $region $file 2>/dev/null` ;
	chomp($data);	
	my ($seq,$from,$to,$name,$coverage) = split(/\t/ , $data);


	my $entry = "<dt>$sample</dt><dd><samp>$coverage</samp></dd>" ;
	chomp($entry);
	printf "    $entry\n";

}

closedir(DIR);

printf "  </dl>\n";

exit 0;
