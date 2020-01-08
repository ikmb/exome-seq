#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Excel::Writer::XLSX;

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

die "Must specify an outfile (--outfile)" unless (defined $outfile);

# Inintiate the XLS workbook
my $workbook = Excel::Writer::XLSX->new($outfile);

# Add a new sheet
my $worksheet = $workbook->add_worksheet();

my $row = 0;

printf "Gene/Transcript/Exon\tMean coverage"

open (my $IN, '<', $infile) or die "FATAL: Can't open file: $infile for reading.\n$!\n";

chomp(my @lines = <$IN>);

# Write a header for this sheet
my @header = ("Exon","Coverage");
&write_xlsx($worksheet, $row, @header);
++$row;

foreach my  $line (@lines) {

	chomp($line);
	my ($chrom,$start,$end,$length,$name,$gc,$mean_coverage,$normalized_coverage,$min_normalized_coverage,$max_normalized_coverage,$min_coverage,$max_coverage,$pct_0x,$read_count) = split("\t", $line);

	if ($mean_coverage < $min_cov) {
		unless (exists($kill{$name})) {
			
			my @ele = ( $name, $mean_coverage );

			&write_xlsx($worksheet, $row, @ele);

			++$row;
		}
	}
	
}

close($IN);

sub write_xlsx{
    my ($worksheet, $tem_row, @ele) = @_;
    for(my $i = 0; $i < @ele; ++$i){
        $worksheet->write( $tem_row, $i, $ele[$i]);
    }
}
