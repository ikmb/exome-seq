include { Strelka ; Strelka_JointCalling } from './../../modules/strelka/main.nf'
include { merge_vcf; vcf_get_sample } from './../../modules/vcf/main.nf'

workflow STRELKA_VARIANT_CALLING {

	take:
	bam
	fasta
	bed
	sample_names

	main:

	if (params.joint_calling) {	
		Strelka_JointCalling(
			bam,
			fasta.collect(),
			bed.collect()
		)
		vcf_merged = Strelka_JointCalling.out.vcf
		vcf_get_sample(
			Strelka_JointCalling.out.vcf.collect(),
			sample_names
		).set { vcf }

	} else {
		Strelka(
			bam,
			fasta.collect(),
			bed.collect()
		)
		vcf = Strelka.out.vcf
		vcf_merged = merge_vcf(
			Strelka.out.vcf.map { m,v ->
				tuple(meta.caller,v)
			}.groupTuple().collect()
		)
	}

	emit:
	vcf = Strelka.out.vcf
	vcf_multi = vcf_merged
}
