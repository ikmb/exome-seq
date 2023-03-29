include { STRELKA_JOINT_CALLING } from './../../modules/strelka/strelka_joint_calling.nf'
include { WHATSHAP } from "./../../modules/whatshap/main.nf"
include { PICARD_SORTVCF  as VCF_GATK_SORT } from "./../../modules/picard/sortvcf"
include { BCFTOOLS_ANNOTATE  as VCF_ADD_HEADER } from "./../../modules/bcftools/annotate"
include { TABIX as VCF_INDEX } from "./../../modules/htslib/tabix"
include { BCFTOOLS_ANNOTATE_DBSNP as VCF_ADD_DBSNP } from "./../../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_VIEW as VCF_FILTER_PASS } from "./../../modules/bcftools/view"
include { BCFTOOLS_MERGE as MERGE_VCF } from "./../../modules/bcftools/merge"
include { GATK_SELECTVARIANTS as VCF_GET_SAMPLE } from "./../../modules/gatk/selectvariants"
include { BCFTOOLS_NORMALIZE as VCF_INDEL_NORMALIZE } from './../../modules/bcftools/normalize'

ch_versions = Channel.from([])

workflow STRELKA_MULTI_CALLING {

	take:
	bams
	bed
	metas
	fasta

	main:

	ch_merged_vcf = Channel.empty()
	ch_phased_multi = Channel.empty()
        ch_vcf = Channel.empty()
    
        bams.map { b,i -> [ [id: "merge"],b,i ] }
        .groupTuple()
        .set { ch_bams }

	STRELKA_JOINT_CALLING(
            ch_bams.map { m,b,i -> [ b,i]},
            bed.collect(),
	    fasta.collect()
        )

	ch_versions = ch_versions.mix(STRELKA_JOINT_CALLING.out.versions)

	VCF_FILTER_PASS(
		STRELKA_JOINT_CALLING.out.vcf.map { v,t -> [ [id: "all", sample_id: "STRELKA_JOINT_CALLING", patient_id: "MergedCallset", variantcaller: "STRELKA"],v,t]}
	)

	ch_versions = ch_versions.mix(VCF_FILTER_PASS.out.versions)

	VCF_GATK_SORT(
		VCF_FILTER_PASS.out.vcf
	)

        ch_merged_vcf = ch_merged_vcf.mix(
		VCF_GATK_SORT.out.vcf
	)

	// Phase Multi-VCF with all samples
	WHATSHAP(
		ch_merged_vcf,
		bams.collect(),
		fasta.collect()
	)
	ch_versions = ch_versions.mix(WHATSHAP.out.versions)

	ch_phased_multi = ch_phased_multi.mix(WHATSHAP.out.vcf)
	
        VCF_GET_SAMPLE(
                ch_merged_vcf.collect(),
		metas
        )

	ch_versions = ch_versions.mix(VCF_GET_SAMPLE.out.versions)

	VCF_INDEL_NORMALIZE(
		VCF_GET_SAMPLE.out.vcf
	)

	ch_versions = ch_versions.mix(VCF_INDEL_NORMALIZE.out.versions)

        VCF_INDEL_NORMALIZE.out.vcf.map { m,v,t ->
                new_meta = m.clone()
                new_meta.variantcaller = "STRELKA"
                tuple(new_meta,v,t)
        }.set { ch_vcf_header }

        VCF_ADD_HEADER(ch_vcf_header)

        single_vcf = ch_vcf.mix(VCF_ADD_HEADER.out.vcf)

	emit:
	versions = ch_versions
	vcf = single_vcf
	vcf_multi = ch_merged_vcf
	vcf_phased_multi = ch_phased_multi
}

