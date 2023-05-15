include { STRELKA_SINGLE_SAMPLE as STRELKA } from './../../modules/strelka/strelka_single.nf'
include { PICARD_SORTVCF  as VCF_GATK_SORT } from "./../../modules/picard/sortvcf"
include { BCFTOOLS_ANNOTATE  as VCF_ADD_HEADER } from "./../../modules/bcftools/annotate"
include { TABIX as VCF_INDEX } from "./../../modules/htslib/tabix"
include { BCFTOOLS_ANNOTATE_DBSNP as VCF_ADD_DBSNP } from "./../../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_VIEW as VCF_FILTER_PASS } from "./../../modules/bcftools/view"
include { BCFTOOLS_MERGE as MERGE_VCF } from "./../../modules/bcftools/merge"
include { GATK_SELECTVARIANTS as VCF_GET_SAMPLE } from "./../../modules/gatk/selectvariants"
include { BCFTOOLS_NORMALIZE as VCF_INDEL_NORMALIZE } from './../../modules/bcftools/normalize'

ch_versions = Channel.from([])
ch_phased_vcf = Channel.from([])
ch_merged_vcf = Channel.empty()
ch_vcf = Channel.empty()

workflow STRELKA_SINGLE_CALLING {

    take:
    bam
    bed
    sample_names
    fasta
    dbsnp

    main:

    STRELKA(
        bam,
        bed,
        fasta
    )

    ch_versions = ch_versions.mix(STRELKA.out.versions)
    
    VCF_INDEX(
        STRELKA.out.vcf
    )
    VCF_FILTER_PASS(
        VCF_INDEX.out.vcf.map { m,v,t ->
            [[
                patient_id: m.patient_id,
                sample_id: m.sample_id,
                variantcaller: "STRELKA"
            ],v,t]
        }
    )

    ch_versions = ch_versions.mix(VCF_FILTER_PASS.out.versions)

        VCF_ADD_DBSNP(
        VCF_FILTER_PASS.out.vcf,
        dbsnp
    )

    ch_versions = ch_versions.mix(VCF_ADD_DBSNP.out.versions)

    VCF_ADD_HEADER(
        VCF_ADD_DBSNP.out.vcf
    )
    single_vcf = ch_vcf.mix(VCF_ADD_HEADER.out.vcf)

    single_vcf.map { m,v,t ->
        new_key = m.sample_id
        tuple(new_key,m,v,t)
    }.set { ch_vcfs_key }

    // Fiddly work-around to determine whether we have 1 or multiple vcfs. No merging when n=1
    VCF_FILTER_PASS.out.vcf.map { m,v,t ->
            [[
                id: "all", 
                sample_id: "UNDEFINED", 
                patient_id: "UNDEFINED", 
                variantcaller: "STRELKA" 
            ],v,t]
    }
    .groupTuple()
    .branch { m,v,t ->
                        single: v.size() == 1
                        multi: v.size() > 1
        }.set { ch_grouped_vcfs }

    MERGE_VCF(
                ch_grouped_vcfs.multi.map { m,v,t -> 
            [ [ 
                id: "all", 
                sample_id: "BCFTOOLS", 
                patient_id: "MERGED_CALLSET", 
                variantcaller: "STRELKA"
            ],v,t] 
        }
    )
    
    ch_versions = ch_versions.mix(MERGE_VCF.out.versions)
        ch_merged_vcf = ch_merged_vcf.mix(MERGE_VCF.out.vcf)

    emit:
    versions = ch_versions
    vcf = single_vcf
    vcf_multi = ch_merged_vcf
    ch_phased_multi = Channel.empty()
}
