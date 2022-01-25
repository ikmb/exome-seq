include { expansion_hunter; expansions2xls } from "./../../modules/expansions/main.nf" params(params)

workflow expansions {

	take:
		bam
		catalog	
	main:
		expansion_hunter(bam,catalog.collect())
		expansions2xls(expansion_hunter.out[0])

	emit:
		expansions = expansions2xls.out
	
}
