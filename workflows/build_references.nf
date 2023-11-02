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
include { VEP_INSTALL_PLUGINS } from "./../modules/vep/install_plugins"
include { DRAGMAP_INSTALL_REFERENCE } from "./../modules/dragmap/install_reference"

fasta_file = file(params.genomes[params.assembly].fasta_ref)

if (params.genomes[params.assembly].report_ref) { report_file = file(params.genomes[params.assembly].report_ref) } else {  report_file = file("${baseDir}/README_PANEL_CALLS.txt") }
if (params.genomes[params.assembly].dragmap_ref) { dragmap_ref = file(params.genomes[params.assembly].dragmap_ref) } else { dragmap_ref = file("${baseDir}/README_PANEL_CALLS.txt") }

if (assembly = "CHM13v2") { version = "vt2t" } else { version = "v38" }

if (params.genomes["refs"][version].dbsnp_ref) { dbsnp_file = params.genomes["refs"][version].dbsnp_ref } else { dbsnp_file = null }
if (params.genomes["refs"][version].mills_ref) { mills_file = params.genomes["refs"][version].mills_ref } else { mills_file = null }
if (params.genomes["refs"][version].g1k_ref) { g1k_file = params.genomes["refs"][version].g1k_ref } else { g1k_file = null }
if (params.genomes["refs"][version].axiom_ref) { axiom_file = params.genomes["refs"][version].axiom_ref } else { axiom_file = null }
if (params.genomes["refs"][version].gnomad_ref) { gnomad_file = params.genomes["refs"][version].gnomad_ref } else { gnomad_file = null }
if (params.genomes["refs"][version].omni_ref) { omni_file = params.genomes["refs"][version].omni_ref } else { omni_file = null }
if (params.genomes["refs"][version].hapmap_ref) { hapmap_file = params.genomes["refs"][version].hapmap_ref } else { hapmap_file = null }

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
).filter { it -> it != null }

ch_report = Channel.fromPath(report_file)

if (params.assembly.contains("alt")) { mode = "ALT" } else { mode = "FULL" }
if (params.genomes[params.assembly].dragmap_ref) { dragmap_build = false } else { dragmap_build = true }
if (fasta_file.getBaseName().contains("GCF_000")) { ncbi = true } else { ncbi = false }

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

ch_fasta = Channel.empty()

workflow BUILD_REFERENCES {

    if ('vep' in tools) {

        // we use the path, not a file, to run the download inside the process
        vep_cache = params.genomes["refs"].vep_ref

        VEP_INSTALL_CACHE(
            vep_cache
        )

        vep_plugins = file(params.genomes["refs"].vep_plugin_ref)
        
        VEP_INSTALL_PLUGINS(
            vep_plugins
        )

    }
    
    // Assembly is gzipped, must unpack first. 
    if (fasta_file.getName().contains(".gz")) {  
        GUNZIP(
            create_genome_channel(fasta_file)
	)

        ch_fasta_gunzip = GUNZIP.out.decompressed

    } else {

        ch_fasta_gunzip = create_genome_channel(fasta_file)

    }

    // assembly is in native NCBI nomenclature, must rename to UCSC convention
    if (ncbi == true) {

        NCBI_ASSEMBLY_CLEAN(
            ch_fasta_gunzip.combine(ch_report)
        )

        ch_fasta = NCBI_ASSEMBLY_CLEAN.out.fasta

    } else {

        RENAME_ASSEMBLY(ch_fasta_gunzip)
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

    // Build dragen index from the assembly fasta file; or download existing index from Illumina
    if (dragmap_build) {
        DRAGMAP_INDEX(
            ch_fasta
        )
    } else {
        DRAGMAP_INSTALL_REFERENCE(
            [[
                assembly: params.assembly
            ],dragmap_ref ]
        )
    }

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

    if (vcf == null) {

    } else {
        def meta = [:]
        meta.assembly     = params.assembly
        meta.id           = file(vcf).getSimpleName()
        meta.file_name    = file(vcf).getName()
        meta.patient_id   = params.assembly
        meta.sample_id    = file(vcf).getSimpleName()

        def array = [ meta, file(vcf) ]
 
        return array
    }
}

