#!/bin/env perl

use Bio::FeatureIO;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

my $usage = qq{
perl ensembl_get_utr.pl
  Getting help:
    [--help]

  Input data
	[--list filename]
	The list of gene names to get BED coordinates for

	[--assembly name]
	Name of the genome assembly to use (GRCH37, hg19, GRCh38)	

  Ouput:
    [--output_file filename]
        The name of the output file. By default the output is the
        standard output
};

my $species = "human";
my $output_file = undef;
my $list = undef;
my $assembly = "GRCh38";
my $help;

GetOptions(
    "help" => \$help,
    "assembly=s" => \$assembly,
    "list=s" => \$list,
    "output_file=s" => \$output_file);
    

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $options = 3306;
my $prefix = "chr";

if ($assembly eq "GRCh37" or $assembly eq "hg19") {
	$options = 3337;
} elsif ($assembly eq "GRCh38") {
	# do nothing
} else {
	exit 1, "Unknown assembly version provided should be one of: hg19, GRCh37 or GRCh38 (default).\n";
}
if ($assembly eq "GRCh37") {
	$prefix = "";
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => $options);

my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );

my $fh = IO::File->new();
$fh->open( $list);

foreach $line (<$fh>) {

	chomp($line);
	
	# Theoretically, one HGNC can map to multiple Genes
	my $genes = $gene_adaptor->fetch_all_by_external_name($line);
	
	my $skip_this = 0;
	
	die "Gene not found ($line)\n" if (scalar @$genes == 0) ;

	foreach my $gene(@$genes) {
		
		# Check if this version of the gene is the reference version
		my $is_reference = $gene->is_reference;

		next if ($is_reference == 0 || $skip_this == 1 || $gene->stable_id =~ /LRG.*/) ;
	
		my $transcript = $gene->canonical_transcript;
		my @exons = @{ $transcript->get_all_translateable_Exons() } ;
		foreach my $exon (@exons) {
			next if (!$exon->is_coding) ;
			my $ref_start = $exon->coding_region_start($transcript);
			my $ref_end = $exon->coding_region_end($transcript);
			if ($ref_start > $ref_end) {
				($ref_start,$ref_end) = ($ref_end,$ref_start);
			}
			my $strand = $exon->strand == 1 ? "+" : "-" ;
			printf $prefix . $gene->seq_region_name . "\t" . $ref_start . "\t" . $ref_end . "\t" . $line . "." . $exon->rank($transcript) . "\t" . 100 . "\t" . $strand . "\n";
		}
					
		# Only need the first occurence of a reference gene
		$skip_this = 1;
	}
}
close ($fh);


