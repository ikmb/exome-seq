include { EXPANSIONS2XLSX } from "./../modules/expansions/expansions2xls"
include { EXPANSION_HUNTER } from  "./../modules/expansions/hunter"

ch_versions = Channel.from([])

workflow EXPANSIONS {

	take:
		bam
		catalog	
	main:
		EXPANSION_HUNTER(bam,catalog.collect())

		ch_versions = ch_versions.mix(EXPANSION_HUNTER.out.versions)

		EXPANSIONS2XLSX(EXPANSION_HUNTER.out[0])

	emit:
		expansions = EXPANSIONS2XLSX.out
		versions = ch_versions	
}
