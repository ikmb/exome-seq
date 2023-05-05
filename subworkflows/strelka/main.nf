include { STRELKA ; STRELKA_JOINTCALLING } from './../../modules/strelka/main.nf'
include { VCF_GATK_SORT ; VCF_ADD_HEADER; VCF_INDEX; VCF_ADD_DBSNP; VCF_FILTER_PASS; MERGE_VCF; VCF_GET_SAMPLE } from './../../modules/vcf/main.nf'
include { VCF_INDEL_NORMALIZE } from './../../modules/bcftools/main.nf'
include { WHATSHAP; WHATSHAP_SINGLE } from './../../modules/whatshap/main.nf'

workflow STRELKA_VARIANT_CALLING {

        ch_merged_vcf = Channel.empty()
        ch_vcf = Channel.empty()

	take:
	bam
	bed
	sample_names

	main:

	STRELKA(
		bam,
		bed.collect()
	)
	VCF_INDEX(STRELKA.out.vcf)
        VCF_FILTER_PASS(VCF_INDEX.out.vcf)
        VCF_ADD_DBSNP(VCF_FILTER_PASS.out.vcf)
	VCF_ADD_HEADER(VCF_ADD_DBSNP.out.vcf.map { meta,v,t ->
			new_meta = [ id: meta.id, sample_id: meta.sample_id, patient_id: meta.patient_id, variantcaller: "STRELKA" ]
			tuple(new_meta,v,t)
		}
	)
	single_vcf = ch_vcf.mix(VCF_ADD_HEADER.out.vcf)

	bam.map { m,b,i -> 
		new_key = m.sample_id
		tuple(new_key,b,i)
	}.set { ch_bams_key }

	single_vcf.map { m,v,t ->
		new_key = m.sample_id
		tuple(new_key,m,v,t)
	}.set { ch_vcfs_key }

	//WHATSHAP_SINGLE(
	//	ch_vcfs_key.join(ch_bams_key).map { n,m,v,t,b,i -> [ m,v,t,b,i ] }
	//)

	// Fiddly work-around to determine whether we have 1 or multiple vcfs. No merging when n=1
	VCF_FILTER_PASS.out.vcf.map { m,v,t ->
		def new_meta = [ id: "all", sample_id: "UNDEFINED", patient_id: "UNDEFINED", variantcaller: "STRELKA" ]
		tuple(new_meta,v,t)
	}
	.groupTuple()
	.branch { m,v,t ->
                        single: v.size() == 1
                        multi: v.size() > 1
        }.set { ch_grouped_vcfs }

	MERGE_VCF(
                ch_grouped_vcfs.multi.map { m,v,t -> [ [ id: "all", sample_id: "Bcftools", patient_id: "MERGED_CALLSET", variantcaller: "STRELKA"],v,t] }
	)
        ch_merged_vcf = ch_merged_vcf.mix(MERGE_VCF.out.vcf)

	emit:
	vcf = single_vcf
	vcf_multi = ch_merged_vcf
	ch_phased_multi = Channel.empty()
}

workflow STRELKA_MULTI_CALLING {

	take:
	bams
	bed
	metas

	main:

	ch_merged_vcf = Channel.empty()
	ch_phased_multi = Channel.empty()
        ch_vcf = Channel.empty()
    
        bams.map { b,i -> [ [id: "merge"],b,i ] }
        .groupTuple()
        .set { ch_bams }

	STRELKA_JOINTCALLING(
            ch_bams.map { m,b,i -> [ b,i]},
            bed.collect()
        )
	VCF_FILTER_PASS(
		STRELKA_JOINTCALLING.out.vcf.map { v,t -> [ [id: "all", sample_id: "STRELKA_JOINT_CALLING", patient_id: "MERGED_CALLSET", variantcaller: "STRELKA"],v,t]}
	)
	VCF_GATK_SORT(
		VCF_FILTER_PASS.out.vcf
	)

        ch_merged_vcf = ch_merged_vcf.mix(
		VCF_GATK_SORT.out.vcf
	)

	// Phase Multi-VCF with all samples
	WHATSHAP(
		ch_merged_vcf,
		bams.collect()
	)
	ch_phased_multi = ch_phased_multi.mix(WHATSHAP.out.vcf)
	
        VCF_GET_SAMPLE(
                ch_merged_vcf.collect(),
		metas
        ).set { svcf }

	VCF_INDEL_NORMALIZE(
		svcf
	)

        VCF_INDEL_NORMALIZE.out.vcf.map { m,v,t ->
                new_meta = m.clone()
                new_meta.variantcaller = "STRELKA"
                tuple(new_meta,v,t)
        }.set { ch_vcf_header }

        VCF_ADD_HEADER(ch_vcf_header)

        single_vcf = ch_vcf.mix(VCF_ADD_HEADER.out.vcf)

	emit:
	vcf = single_vcf
	vcf_multi = ch_merged_vcf
	vcf_phased_multi = ch_phased_multi
}

