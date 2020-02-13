#!/usr/bin/env perl

use Cwd;
use Getopt::Long;

my $usage = qq{
  Getting help:
    [--help] 
    [--region]
	A BED file with the target region to check
    [--fasta]
	The genome sequence in FASTA format to decode CRAM compression (optional)

  Ouput:
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $region = undef;
my $fasta = undef;
my $help;

GetOptions(
    "help" => \$help,
    "fasta=s" => \$fasta,
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
id: 'Coverage of the SRY gene'
section_name: 'Sex check'
plot_type: 'html'
description: ' shows the coverage of the SRY gene to determine sex of the sample. Males will have a high coverage.'
data: |\n
  <dl class="dl-horizontal">
);

printf $header . "\n";

opendir(DIR, $dir) or die $!;

while (my $file = readdir(DIR) ) {

	next unless ($file =~ m/\.bam$/ || $file =~ m/\.cram$/ );
	
	my $sample = (split /\./, $file )[0] ;

	my $command = "";

	if (defined $fasta) {
		$command = "samtools bedcov --reference $fasta $region $file 2>/dev/null" ;
	} else {
		$command = "samtools bedcov $region $file 2>/dev/null" ;
	}

	my $data = `$command` ;
	chomp($data);	
	my ($seq,$from,$to,$name,$coverage) = split(/\t/ , $data);

	# The length of SRY is 827bp, so we divide the total coverage base count of the bedcov tool by that

	my $normalized_coverage = $coverage/827 ;

	my $entry = "<dt>$sample</dt><dd><samp>$normalized_coverage</samp></dd>" ;
	chomp($entry);
	printf "    $entry\n";

}

closedir(DIR);

printf "  </dl>\n";

exit 0;
