include { GATK_HAPLOTYPECALLER as GATK_HAPLOTYPECALLER_SINGLE; GATK_HAPLOTYPECALLER as GATK_HAPLOTYPECALLER_GVCF} from './../modules/gatk/haplotypecaller'
include { GATK_COMBINEGVCFS } from './../modules/gatk/combinegvcfs'
include { GATK_GENOTYPEGVCFS } from './../modules/gatk/genotypegvcfs'
include { GATK_GENOMICSDBIMPORT } from './../modules/gatk/genomicsdbimport'
include { GATK_MAKESITESONLYVCF } from './../modules/gatk/sitesonly'
include { GATK_VARIANTFILTRATION } from './../modules/gatk/variantfiltration'
include { GATK_VARIANTRECALIBRATOR as GATK_INDEL_RECALIBRATOR; GATK_VARIANTRECALIBRATOR as GATK_SNP_RECALIBRATOR } from './../modules/gatk/variantrecalibrator'
include { GATK_APPLYVQSR as GATK_INDEL_VQSR; GATK_APPLYVQSR as GATK_SNP_VQSR } from './../modules/gatk/applyvqsr'
include { GATK_CNNSCOREVARIANTS } from './../modules/gatk/cnnscorevariants'
include { GATK_MERGEVCFS } from './../modules/gatk/mergevcfs'
include { GATK_FILTERVARIANTTRANCHES } from './../modules/gatk/filtervarianttranches'

ch_vcf_multi = Channel.from([])
ch_vcf_single = Channel.from([])

workflow GATK_VARIANT_CALLING {

	take:
	bam
	intervals
	fasta
	known_snps
	known_snps_tbi
	known_indels
	known_indels_tbi

	main:
	
	if (params.joint_calling) {

		// Produce gVCFs
		GATK_HAPLOTYPECALLER_GVCF(
			bam,
			intervals.collect(),
			"gvcf",
			fasta.collect()
		)

		// Combine all gVCFs into one multi-sample gVCF
		GATK_COMBINEGVCFS(
			GATK_HAPLOTYPECALLER_GVCF.out.vcf.map { m,v,t -> v }.collect(),
				GATK_HAPLOTYPECALLER_GVCF.out.vcf.map { m,v,t -> t }.collect(),
			intervals.collect(),
			fasta.collect()
        )

		// Call genotypes from gVCF(s)
		GATK_GENOTYPEGVCFS(
			GATK_COMBINEGVCFS.out.gvcf,
			intervals.collect(),
			fasta.collect()
        )

		GATK_GENOTYPEGVCFS.out.vcf.map { v,t ->
                        new_meta = [ id: "all", sample_id: "UNDEFINED", patient_id: "UNDEFINED", variantcaller: "GATK" ]
                        tuple(new_meta,v,t)
                }.set { ch_variants_pass }

		// Hard-filter variants on excess heterosygosity
		GATK_VARIANTFILTRATION(
			ch_variants_pass
		)

		ch_multi_vcf_filtered = GATK_VARIANTFILTRATION.out.vcf.map { v,t ->
			def new_meta = [ id: "all", sample_id: "GATK", patient_id: "MergedCallset", variantcaller: "GATK"]
			tuple(new_meta,v,t)
		}
		
		// Produce a sites-only vcf to speed up variant recalibration
		GATK_MAKESITESONLYVCF(
			ch_multi_vcf_filtered
		)

		// Compute indel recalibration
		GATK_INDEL_RECALIBRATOR(
			GATK_MAKESITESONLYVCF.out.vcf,
			"INDEL"
		)	

		// Compute snp recalibration
		GATK_SNP_RECALIBRATOR(
			GATK_MAKESITESONLYVCF.out.vcf,
			"SNP"
		)

		// Apply indel recalibration
		GATK_INDEL_VQSR(
			ch_multi_vcf_filtered,
			GATK_INDEL_RECALIBRATOR.out.recal,
			GATK_INDEL_RECALIBRATOR.out.tranches,
			"INDEL"
		)

		// Apply snp recalibration
		GATK_SNP_VQSR(
			ch_multi_vcf_filtered,
			GATK_SNP_RECALIBRATOR.out.recal,
			GATK_SNP_RECALIBRATOR.out.tranches,
			"SNP"
		)

		// Merge indel snd snp recalibrated vcfs
		GATK_MERGEVCFS(
			GATK_SNP_VQSR.out.vcf.map { m,v,t -> [ m,v ] }.join(
				GATK_INDEL_VQSR.out.vcf.map { m,v,t -> [ m,v ] }
			),
			GATK_SNP_VQSR.out.vcf.map { m,v,t -> [ m,t ] }.join(
				GATK_INDEL_VQSR.out.vcf.map { m,v,t -> [ m,t ] }
			)
		)

		ch_vcf_multi = ch_vcf_multi.mix(GATK_MERGEVCFS.out.vcf)
	
	} else {

		// Call single samples
		GATK_HAPLOTYPECALLER_SINGLE(
			bam,
			intervals.collect(),
			"single",
			fasta.collect()
		)

		// Filter VCF using GATK neural networks
		GATK_CNNSCOREVARIANTS(
			GATK_HAPLOTYPECALLER_SINGLE.out.vcf.join(
				GATK_HAPLOTYPECALLER_SINGLE.out.bam
			),
			intervals.collect(),
			fasta.collect()
		)

		GATK_FILTERVARIANTTRANCHES(
			GATK_CNNSCOREVARIANTS.out.vcf,
			known_snps.collect(),
			known_snps_tbi.collect(),
			known_indels.collect(),
			known_indels_tbi.collect()
		)
		ch_vcf_single = ch_vcf_single.mix(GATK_FILTERVARIANTTRANCHES.out.vcf)

	}

	emit:
	vcf = ch_vcf_single
	vcf_multi = ch_vcf_multi

	//vcf = VCF_GET_SAMPLE.out.vcf
	//vcf_multi = GATK_GENOTYPEGVCFS.out.vcf
}
