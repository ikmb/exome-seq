include { NCBI_ASSEMBLY_CLEAN } from "./../modules/helper/ncbi_assembly_clean"
include { BWA2_MEM_INDEX } from "./../modules/bwa2/index"
include { BWA_MEM_INDEX } from "./../modules/bwa/index"
include { DRAGMAP_INDEX } from "./../modules/dragmap/index"
include { GUNZIP } from "./../modules/gunzip" 
include { SAMTOOLS_DICT } from "./../modules/samtools/dict"
include { SAMTOOLS_FAIDX } from "./../modules/samtools/faidx"
include { BGZIP_INDEX } from "./../modules/htslib/bgzip_index"
include { TABIX } from "./../modules/htslib/tabix"
include { BCFTOOLS_RENAME_CHROMOSOMES } from "./../modules/bcftools/rename_chromosomes"
include { PICARD_UPDATE_SEQUENCE_DICTIONARY } from "./../modules/picard/update_sequence_dictionary"
include { BCFTOOLS_REHEADER } from "./../modules/bcftools/reheader"
include { RENAME_ASSEMBLY } from "./../modules/helper/rename_assembly"
include { VEP_INSTALL_CACHE } from "./../modules/vep/install_cache"

fasta_file = file(params.genomes[params.assembly].fasta_ref)

if (params.genomes[params.assembly].report_ref) { report_file = file(params.genomes[params.assembly].report_ref) } else {  report_file = file("${baseDir}/README_PANEL_CALLS.txt") }

dbsnp_file = file(params.genomes["refs"].dbsnp_ref)
mills_file = file(params.genomes["refs"].mills_ref)
g1k_file = file(params.genomes["refs"].g1k_ref)
axiom_file = file(params.genomes["refs"].axiom_ref)
gnomad_file = file(params.genomes["refs"].gnomad_ref)
omni_file = file(params.genomes["refs"].omni_ref)
hapmap_file = file(params.genomes["refs"].hapmap_ref)

ch_variants = Channel.fromList(
	[ 
		create_vcf_channel(mills_file),
		create_vcf_channel(g1k_file),
		create_vcf_channel(axiom_file),
		create_vcf_channel(gnomad_file),
		create_vcf_channel(omni_file),
		create_vcf_channel(dbsnp_file),
        create_vcf_channel(hapmap_file)
	]
)

ch_report = Channel.fromPath(report_file)

if (params.assembly.contains("alt")) {
    mode = "ALT"
} else {
    mode = "FULL"
}

if (fasta_file.getBaseName().contains("GCF_000")) {
    ncbi = true
} else {
    ncbi = false
}

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

ch_fasta = Channel.empty()

workflow BUILD_REFERENCES {

    if ('vep' in tools) {

        vep_cache = file(params.genomes["refs"].vep_ref)

        VEP_INSTALL_CACHE(
            vep_cache
        )

        vep_plugins = file(params.genomes["refs"].vep_plugin_ref)
        
        VEP_INSTALL_PLUGINS(
            vep_plugins
        )

    }
    
    GUNZIP(
        create_genome_channel(fasta_file)
    )

    if (ncbi == true) {

        NCBI_ASSEMBLY_CLEAN(
            GUNZIP.out.decompressed.combine(ch_report)
        )

        ch_fasta = NCBI_ASSEMBLY_CLEAN.out.fasta

    } else {

        RENAME_ASSEMBLY(GUNZIP.out.decompressed)
        ch_fasta = RENAME_ASSEMBLY.out.fasta

    }  

    SAMTOOLS_DICT(
        ch_fasta
    )

    SAMTOOLS_FAIDX(
        ch_fasta
    )

    BWA2_MEM_INDEX(
        ch_fasta
    )

    BWA_MEM_INDEX(
        ch_fasta
    )

    DRAGMAP_INDEX(
        ch_fasta
    )

    TABIX(
        ch_variants
    )

    BCFTOOLS_REHEADER(
        TABIX.out.vcf,
        SAMTOOLS_FAIDX.out.fai.collect()
    )

    PICARD_UPDATE_SEQUENCE_DICTIONARY(
        BCFTOOLS_REHEADER.out.vcf,
        SAMTOOLS_DICT.out.dict.collect()
    )

}

def create_genome_channel(genome) {
    def meta = [:]
    meta.assembly     = params.assembly
    meta.id           = file(genome).getSimpleName()

    def array = [ meta, genome ]

    return array
}

def create_vcf_channel(vcf) {
    def meta = [:]
    meta.assembly     = params.assembly
    meta.id           = file(vcf).getSimpleName()
    meta.file_name    = file(vcf).getName()
    meta.patient_id   = params.assembly
    meta.sample_id    = file(vcf).getSimpleName()

    def array = [ meta, vcf ]

    return array
}
