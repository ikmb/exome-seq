#!/bin/env perl

use Bio::FeatureIO;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

my $usage = qq{
perl ensembl_get_utr.pl
  Getting help:
    [--help]

  For the species
  	[--species]
  	Names of species to use
	
  Ouput:
    [--output_file filename]
        The name of the output file. By default the output is the
        standard output
};

my $species = undef;
my $output_file = undef;
my $list = undef;
my $version = "GRCh38";
my $help;

GetOptions(
    "help" => \$help,
    "species=s" => \$species,
	"version" => \$version,
    "list=s" => \$list,
    "output_file=s" => \$output_file);
    

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $options = "";
my $prefix = "chr";

if ($version == "GRCh37" || $version == "hg19") {
	$options = "-port 3337";
} elsif ($version == "GRCh38") {
	# do nothing
} else {
	exit 1, "Unknown assembly version provided should be one of: hg19, GRCh37 or GRCh38 (default).\n";
}
if ($version == "GRCh37") {
	$prefix = "";
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
	$options);

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
		my @exons = @{ $transcript->get_all_translatable_Exons() } ;
		foreach my $exon (@exons) {
			my $strand = $exon->strand == 1 ? "+" : "-" ;
			printf $prefix . $gene->seq_region_name . "\t" . $exon->start . "\t" . $exon->end . "\t" . $line . "." . $exon->rank($transcript) . "\t" . 100 . "\t" . $strand . "\n";
		}
					
		# Only need the first occurence of a reference gene
		$skip_this = 1;
	}
}
close ($fh);


