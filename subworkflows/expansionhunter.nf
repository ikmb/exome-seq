include { EXPANSIONS2XLSX } from "./../modules/expansions/expansions2xls"
include { EXPANSION_HUNTER } from  "./../modules/expansions/hunter"

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
