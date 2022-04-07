include { cnvkit_ref_decompress; cnvkit_autobin ; cnvkit_ref_to_targets ; cnvkit_coverage ;  cnvkit_process ; cnvkit_segmetrics ; cnvkit_call ; cnvkit_genemetrics ;  cnvkit_breaks ; cnvkit_export ; cnvkit_plots } from "./../../modules/cnvkit/main.nf" params(params)


workflow CNVKIT {

	take:
		bed
		bam
		cnn_gz
	main:
		cnvkit_ref_decompress(cnn_gz)
		cnvkit_ref_to_targets(cnvkit_ref_decompress.out)
		cnvkit_coverage(bam,cnvkit_ref_to_targets.out.collect())
		cnvkit_process(cnvkit_coverage.out[1],cnvkit_ref_decompress.out.collect())
		cnvkit_segmetrics(cnvkit_process.out)
		cnvkit_call(cnvkit_segmetrics.out)
		cnvkit_genemetrics(cnvkit_call.out)
		cnvkit_breaks(cnvkit_call.out)
		cnvkit_export(cnvkit_call.out)
		cnvkit_plots(cnvkit_call.out)

	emit:
		bed = cnvkit_export.out[0]
		vcf = cnvkit_export.out[1]				

}
