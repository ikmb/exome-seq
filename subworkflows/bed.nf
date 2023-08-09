include { INTERVAL_TO_BED } from "./../modules/picard/intervallisttools" 
include { BEDTOOLS_SLOP } from "./../modules/bedtools/slop"
include { BGZIP_INDEX as BED_TO_BEDGZ } from "./../modules/htslib/bgzip_index"

workflow CONVERT_BED {
    
    take:
        intervals
        fasta

    main:
        INTERVAL_TO_BED(intervals)
        BEDTOOLS_SLOP(
            INTERVAL_TO_BED.out.bed,
            params.interval_padding,
            fasta
        )
        BED_TO_BEDGZ(BEDTOOLS_SLOP.out.bed)

    emit:
        bed_padded = BEDTOOLS_SLOP.out.bed
        bed = INTERVAL_TO_BED.out.bed
        bed_gz = BED_TO_BEDGZ.out.bedgz

}
