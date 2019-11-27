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
my $min_cov = 30;
my $ban = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "min_cov=i" => \$min_cov,
    "ban=s" => \$ban,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

my %kill;

if (defined $ban) {
	open my $ban_list , '<', $ban;
	chomp(my @lines = <$ban_list>);
	foreach my $line (@lines) {
		$kill{$line} = 1;
	}

	close $ban_list;
}

my $header = qq(
id: 'low coverage targets'\n
section_name: 'ikmb/exome-seq low coverage targets'\n
plot_type: 'html'\n
description: ' lists all defined exome targets that fall below $min_cov X coverage.'\n
data: |\n
  <dl class="dl-horizontal">
);

printf $header . "\n";

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

chomp(my @lines = <$IN>);

my $header = shift @lines;

foreach my  $line (@lines) {

	chomp($line);
	my ($chrom,$start,$end,$length,$name,$gc,$mean_coverage,$normalized_coverage,$min_normalized_coverage,$max_normalized_coverage,$min_coverage,$max_coverage,$pct_0x,$read_count) = split("\t", $line);

	if ($mean_coverage < $min_cov) {
		unless (exists($kill{$name})) {
			my $entry = "<dt>$name</dt><dd><samp>$mean_coverage</samp></dd>" ;
			chomp($entry);
        		printf "    $entry\n";
		}
	}
	
}

close($IN);

printf "  </dl>\n";
