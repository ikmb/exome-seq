process BGZIP_INDEX {

	container 'quay.io/biocontainers/htslib:1.16--h6bc39ce_0'
        executor 'local'

        input:
        path(bed)

        output:
        tuple path(bed_gz),path(bed_gz_tbi)

        script:
        bed_gz = bed.getBaseName() + ".bed.gz"
        bed_gz_tbi = bed_gz + ".tbi"

        """
                bgzip $bed
                tabix -p bed $bed_gz
        """
}


