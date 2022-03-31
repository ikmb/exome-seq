include { EXPANSION_HUNTER; EXPANSIONS2XLSX } from "./../../modules/expansions/main.nf" params(params)

workflow EXPANSIONS {

	take:
		bam
		catalog	
	main:
		EXPANSION_HUNTER(bam,catalog.collect())
		EXPANSIONS2XLSX(EXPANSION_HUNTER.out[0])

	emit:
		expansions = EXPANSIONS2XLSX.out
	
}
