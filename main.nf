#!/usr/bin/env nextflow

/**
IKMB Diagnostic Exome Pipeline

This Pipeline performs one of two workflows to generate variant calls and effect predictions
using either the GATK processing chain or Freebayes
**/

// #############
// INPUT OPTIONS
// #############

inputFile = file(params.samples)

// Available tool chains

valid_tools = [ "freebayes", "gatk3", "gatk4" ]
params.tool = "freebayes"

if (valid_tools.contains(params.tool) == false) {
   exit 1; "Specified an unknown tool chain, please consult the documentation for valid tool chains."
}

// EXOMISER input data
params.hpo = false
params.ped = false
params.omim = false

// This will eventually enable switching between multiple assembly versions
// Currently, only hg19 has all the required reference files available
params.assembly = "hg19_clinical"

if (params.genomes.containsKey(params.assembly) == false) {
   exit 1, "Specified unknown genome assembly, please consult the documentation for valid assemblies."
}

REF = file(params.genomes[ params.assembly ].fasta)
REF_CHR_ROOT = params.genomes[ params.assembly ].per_sequence_root
DICT = file(params.genomes[params.assembly ].dict)

DBSNP = file(params.genomes[ params.assembly ].dbsnp )
G1K = file(params.genomes[ params.assembly ].g1k )
GOLD1 = file(params.genomes[ params.assembly ].gold )
OMNI = file(params.genomes[ params.assembly ].omni )
HAPMAP = file(params.genomes[ params.assembly ].hapmap )
EXAC = file(params.genomes[ params.assembly ].exac )
CADD = file(params.genomes[ params.assembly ].cadd )
ANNOVAR_DB = file(params.genomes[ params.assembly ].annovar_db )

VEP_CACHE = params.vep_cache

// Location of applications used
GATK = file(params.gatk_jar)
PICARD = file(params.picard_jar)
OUTDIR = file(params.outdir)

// Allow for custom freebayes filter options
params.freebayes_options = "--min-alternate-fraction 0.2 --min-base-quality 20 --min-alternate-qsum 90"
freebayes_options = params.freebayes_options

// Available exome kits
if (params.genomes[params.assembly].kits.containsKey(params.kit) == false) {
   exit 1, "Specified unknown Exome kit, please consult the documentation for valid kits."
}

TARGETS= params.genomes[params.assembly].kits[ params.kit ].targets
BAITS= params.genomes[params.assembly].kits[ params.kit ].baits

// Determine valid intervals for parallel processing from the exome target file
chromosomes =  []
file(TARGETS).eachLine { line ->
	elements = line.trim().split("\t")
	seq = elements[0].trim()
	if (seq =~ /^chr.*/) {
		if (chromosomes.contains(seq) == false)	{
			chromosomes << seq
		}
	}
}

// We add 17 reference exome gVCFs to make sure that variant filtration works
// These are in hg19 so need to be updated to other assemblies if multiple assemblies are to be supported

calibration_exomes = file(params.calibration_exomes)
calibration_samples_list = file(params.calibration_exomes_samples)

// GATK4 requires a specific name extension...
if (params.tool == "gatk4") {
	calibration_samples_list = file(params.calibration_exomes_samples_args)
}

calibration_vcfs = [ ]
file(calibration_exomes).eachLine { line ->
        location = line.trim()
        calibration_vcfs << location
}

// Whether to send a notification upon workflow completion
params.email = false

// Whether to use a local scratch disc
use_scratch = params.scratch

// Trimmomatic options
TRIMMOMATIC=file(params.trimmomatic)

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters

logParams(params, "pipeline_parameters.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "IKMB Diagnostic Exome pipeline v${VERSION}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "Tool Chain:			${params.tool}"
log.info "Adapter sequence used:	${adapters}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

// Read sample file 
Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  readPairsTrimmomatic }

process runTrimmomatic {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    //publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/Trimmomatic/", mode: 'copy'

    scratch use_scratch

    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, date, fastqR1, fastqR2 from readPairsTrimmomatic

    output:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, date, center, file("${libraryID}_R1.paired.fastq.gz"),file("${libraryID}_R2.paired.fastq.gz") into inputBwa, readPairsFastQC

    script:

    """
	java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads 8 $fastqR1 $fastqR2 ${libraryID}_R1.paired.fastq.gz ${libraryID}.1U.fastq.gz ${libraryID}_R2.paired.fastq.gz ${libraryID}.2U.fastq.gz ILLUMINACLIP:${TRIMMOMATIC}/adapters/${adapters}:2:30:10:3:TRUE LEADING:${leading} TRAILING:${trailing} SLIDINGWINDOW:${slidingwindow} MINLEN:${minlen}
     """
}

process runBWA {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    //publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/", mode: 'copy'

    scratch use_scratch
	
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa
    
    output:
    set indivID, sampleID, file(outfile) into runBWAOutput
    
    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"	

    """
	bwa mem -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t 16 ${REF} $left $right | samtools sort -O bam - > $outfile
    """	
}

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])

process mergeBamFiles_bySample {

        tag "${indivID}|${sampleID}"
	
	input:
        set indivID, sampleID, file(aligned_bam_list) from runBWAOutput_grouped_by_sample

	output:
	set indivID,sampleID,file(merged_bam) into mergedBamFile_by_Sample

	script:
	merged_bam = sampleID + "merged.bam"

	"""
		java -jar -Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=tmp/ -jar ${PICARD} MergeSamFiles \
			INPUT=${aligned_bam_list.join(' INPUT=')} \
			OUTPUT=${merged_bam} \
			CREATE_INDEX=false \
			CREATE_MD5_FILE=false \
			SORT_ORDER=coordinate
	"""
}

process runMarkDuplicates {

	tag "${indivID}|${sampleID}"
        publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/MarkDuplicates", mode: 'copy'

        // scratch use_scratch

        input:
        set indivID, sampleID, file(merged_bam) from mergedBamFile_by_Sample

        output:
        set indivID, sampleID, file(outfile_bam),file(outfile_bai) into MarkDuplicatesOutput, BamForMultipleMetrics, runPrintReadsOutput_for_OxoG_Metrics, runPrintReadsOutput_for_HC_Metrics, BamForDepthOfCoverage
	file(outfile_bam) into FreebayesBamInput
	file(outfile_bai) into FreebayesBaiInput
	file(outfile_md5) into MarkDuplicatesMD5
	file(outfile_metrics) into DuplicatesOutput_QC

        script:
        outfile_bam = sampleID + ".dedup.bam"
        outfile_bai = sampleID + ".dedup.bai"
	outfile_md5 = sampleID + ".dedup.bam.md5"

        outfile_metrics = sampleID + "_duplicate_metrics.txt"

	"""
        	java -Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=5 -Djava.io.tmpdir=tmp/ -jar ${PICARD} MarkDuplicates \
                	INPUT=${merged_bam} \
	                OUTPUT=${outfile_bam} \
        	        METRICS_FILE=${outfile_metrics} \
                        CREATE_INDEX=true \
                        TMP_DIR=tmp && md5sum ${outfile_bam} > ${outfile_md5}
			
	"""

}

// *******************
// FREEBAYES WORKFLOW
// *******************

if (params.tool == "freebayes") {

	process runFreebayes {

		tag "ALL|${chr}"
                // publishDir "${OUTDIR}/${params.tool}/Variants/ByChromosome", mode: 'copy'

		input:
		file(bam_files) from FreebayesBamInput.collect()
		file(bai_files) from FreebayesBaiInput.collect()
		each chr from chromosomes
	
		output:
		file(vcf) into outputFreebayes
		
		script:
		vcf = "freebayes.${chr}.vcf"

		"""
			freebayes-parallel <(ruby $baseDir/bin/bed2regions $chr < $TARGETS) ${task.cpus} -f ${REF} $freebayes_options ${bam_files} > ${vcf}
		"""
	
	}

	process runConcatVcf {

		tag "ALL"
                // publishDir "${OUTDIR}/${params.tool}/Variants/Freebayes", mode: 'copy'

		input:
		file(vcf_files) from outputFreebayes.collect()
		
		output:
		file(vcf_merged) into outputVcfMerged

		script:
		vcf_merged = "freebayes.merged.vcf"

		"""
			vcf-concat ${vcf_files.join(" ")} | vcf-sort > $vcf_merged
		"""
	}

	process runFilterVcf {

		tag "ALL"
		// publishDir "${OUTDIR}/${params.tool}/Variants/Freebayes", mode: 'copy'	

		input:
		file(vcf) from outputVcfMerged

		output:
		set file(vcf_filtered),file(vcf_filtered_index) into inputLeftNormalize

		script: 
		vcf_filtered = "freebayes.merged.filtered.vcf.gz"
		vcf_filtered_index = vcf_filtered + ".tbi"

		"""
			vcffilter -f "QUAL > 20" ${vcf} | bgzip > ${vcf_filtered}
			tabix ${vcf_filtered}
		"""
	}

// *********************
// GATK4 WORKFLOW
// *********************
} else if (params.tool == "gatk4") {

	// ------------------------------------------------------------------------------------------------------------
	//
	// Perform base quality score recalibration (BQSR) including
	// 1) Generate a recalibration table
	// 2) Generate a new table after applying recalibration
	// 3) Compare differences between recalibration tables
	// 4) Apply recalibration
	//
	// ------------------------------------------------------------------------------------------------------------

	process runBaseRecalibrator {

    		tag "${indivID}|${sampleID}"
    		// publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Processing/BaseRecalibrator/", mode: 'copy'
	    
    		input:
    		set indivID, sampleID, dedup_bam, dedup_bai from MarkDuplicatesOutput
    
    		output:
    		set indivID, sampleID, dedup_bam, file(recal_table) into runBaseRecalibratorOutput
    
    		script:
    		recal_table = sampleID + "_recal_table.txt" 
       
    		"""
			gatk --java-options "-Xmx${task.memory.toGiga()}G" BaseRecalibrator \
			--reference ${REF} \
			--input ${dedup_bam} \
			--known-sites ${GOLD1} \
			--known-sites ${DBSNP} \
        	        --known-sites ${G1K} \
			--output ${recal_table}
		"""
	}

	process runApplyBQSR {

		tag "${indivID}|${sampleID}"
    		publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/", mode: 'copy'

    		scratch use_scratch
	    
    		input:
    		set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput 

    		output:
    		set indivID, sampleID, file(outfile_bam), file(outfile_bai) into runPrintReadsOutput_for_Multiple_Metrics,inputHCSample
    		set indivID, sampleID, realign_bam, recal_table into runPrintReadsOutput_for_PostRecal
		set indivID, outfile_md5 into BamMD5
            
    		script:
    		outfile_bam = sampleID + ".clean.bam"
    		outfile_bai = sampleID + ".clean.bai"
		outfile_md5 = sampleID + ".clean.bam.md5"
           
    		"""
                gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyBQSR \
                --reference ${REF} \
                --input ${realign_bam} \
                -bqsr ${recal_table} \
                --output ${outfile_bam} \
                -OBM true
    		"""
	}    

	process runBaseRecalibratorPostRecal {

		tag "${indivID}|${sampleID}"
    		// publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Processing/BaseRecalibratorPostRecal/", mode: 'copy'
	    
    		input:
    		set indivID, sampleID, realign_bam, recal_table from runPrintReadsOutput_for_PostRecal
    
    		output:
    		set indivID, sampleID, recal_table, file(post_recal_table) into runBaseRecalibratorPostRecalOutput_Analyze
        
    		script:
    		post_recal_table = sampleID + "_post_recal_table.txt" 
   
    		"""
			gatk --java-options "-Xmx25G" BaseRecalibrator \
			--reference ${REF} \
			--input ${realign_bam} \
			--known-sites ${GOLD1} \
			--known-sites ${DBSNP} \
			--output ${post_recal_table}
		"""
	}	

	process runAnalyzeCovariates {

    		tag "${indivID}|${sampleID}"
    		publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Processing/AnalyzeCovariates/", mode: 'copy'
	    
    		input:
    		set indivID, sampleID, recal_table, post_recal_table from runBaseRecalibratorPostRecalOutput_Analyze

		output:
		set indivID, sampleID, recal_plots into runAnalyzeCovariatesOutput
	    
    		script:
    		recal_plots = sampleID + "_recal_plots.pdf" 

    		"""
			gatk --java-options "-Xmx5G" AnalyzeCovariates \
				--before-report-file ${recal_table} \
				--after-report-file ${post_recal_table} \
				--plots-report-file ${recal_plots}
		"""
	}    

	// Call variants on a per-sample basis

	process runHCSample {

  		tag "${indivID}|${sampleID}"
  		// publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Variants/HaplotypeCaller/perChromosome" , mode: 'copy'

  		input: 
  		set indivID,sampleID,file(bam),file(bai) from inputHCSample

  		output:
  		file(vcf) into outputHCSample
                file(vcf_index) into outputHCSampleIndex

  		script:
  
  		vcf = sampleID + ".raw_variants.g.vcf.gz"
		vcf_index = vcf + ".tbi"

  		"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" HaplotypeCaller \
			-R $REF \
			-I ${bam} \
			-L $TARGETS \
			-L chrM \
			--genotyping-mode DISCOVERY \
			--emit-ref-confidence GVCF \
			-OVI true \
    			--output $vcf \
			--native-pair-hmm-threads ${task.cpus} &> log.txt \
  		"""
	}

	// Import individual vcf files into a GenomicsDB database on a per chromosome basis
	// From here on all samples are in the same file
	process runGenomicsDBImport  {

		tag "${chr} - using 17 IKMB reference exomes for calibration"
                // publishDir "${OUTDIR}/${params.tool}/Variants/JoinedGenotypes/PerRegion"

		input:
                file(vcf_list) from outputHCSample.collect()
		file(index_list) from outputHCSampleIndex.collect()

		each chr from chromosomes

		output:
                set chr,file(genodb) into inputJoinedGenotyping

		script:
		genodb = "genodb_${chr}"

		"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" GenomicsDBImport  \
			--variant ${vcf_list.join(" --variant ")} \
			--reference $REF \
			--intervals $chr \
			--genomicsdb-workspace-path $genodb
		"""

	}

	// Perform genotyping on a per chromosome basis

	process runJoinedGenotyping {
  
  		tag "${chr} - using 17 IKMB reference exomes for calibration"
  		// publishDir "${OUTDIR}/${params.tool}/Variants/JoinedGenotypes/PerRegion"
  
  		input:
  		set chr,file(genodb) from inputJoinedGenotyping
  
  		output:
  		file(gvcf) into inputCombineVariantsFromGenotyping
  
  		script:
  
  		gvcf = "genotypes." + chr + ".gvcf.gz"
  
  		"""
	 	gatk --java-options "-Xmx${task.memory.toGiga()}G" GenotypeGVCFs \
			--reference $REF \
			--dbsnp $DBSNP \
			-V gendb://${genodb} \
                	--output $gvcf \
                        -G StandardAnnotation \
			-OVI true
  		"""
	}

	// Merging the scatter-gather VCF files into one file

	process combineVariantsFromGenotyping {
		tag "ALL - using 17 IKMB reference exomes for calibration"
		// publishDir "${OUTDIR}/${params.tool}/Variants/JoinedGenotypes"

		input:
		file(vcf_files) from inputCombineVariantsFromGenotyping.collect()

		output:
		file(gvcf) into (inputRecalSNP , inputRecalIndel)

		script:

		gvcf = "genotypes.merged.vcf.gz"

		"""
        		vcf-concat ${vcf_files.join(" ")} | vcf-sort | bgzip > $gvcf
		"""
        }

	process runRecalibrationModeSNP {

  		tag "ALL"
  		// publishDir "${OUTDIR}/${params.tool}/Variants/Recal"

  		input:
  		file(vcf) from inputRecalSNP

  		output:
	  	set file(recal_file),file(tranches),file(rscript),file(snp_file) into inputRecalSNPApply

  		script:
		snp_file = "genotypes.merged.snps.vcf.gz"
  		recal_file = "genotypes.recal_SNP.recal"
  		tranches = "genotypes.recal_SNP.tranches"
  		rscript = "genotypes.recal_SNP.R"

  		"""

		gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants \
			-R $REF \
			-V $vcf \
			-select-type SNP \
			-O $snp_file

		gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
			-R $REF \
			-V $snp_file \
                	-O $recal_file \
        	        --tranches-file $tranches \
	                --rscript-file $rscript \
			-an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP \
        	        -mode SNP \
			--resource hapmap,known=false,training=true,truth=true,prior=15.0:$HAPMAP \
			--resource omni,known=false,training=true,truth=true,prior=12.0:$OMNI \
			--resource 1000G,known=false,training=true,truth=false,prior=10.0:$G1K \
			--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
  		"""
	}

	process runRecalibrationModeIndel {

  		tag "ALL"
  		// publishDir "${OUTDIR}/${params.tool}/Variants/Recal"

  		input:
  		file(vcf) from inputRecalIndel

  		output:
  		set file(recal_file),file(tranches),file(rscript),file(indel_file) into inputRecalIndelApply

  		script:
		indel_file = "genotypes.merged.indel.vcf.gz"
  		recal_file = "genotypes.recal_Indel.recal"
  		tranches = "genotypes.recal_Indel.tranches"
		rscript = "genotypes.recal_Indel.R"

  		"""
		
		gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants \
			-R $REF \
			-V $vcf \
			-select-type INDEL \
			-select-type MIXED \
			-select-type MNP \
			-select-type SYMBOLIC \
			-select-type NO_VARIATION \
			-O $indel_file

	        gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
        	        -R $REF \
	                -V $indel_file \
                	-O $recal_file \
        	        --tranches-file $tranches \
	                --rscript-file $rscript \
                	-an MQ -an MQRankSum -an SOR -an ReadPosRankSum -an FS  \
        	        -mode INDEL \
	                --resource mills,known=false,training=true,truth=true,prior=15.0:$GOLD1 \
                	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
  		"""

	}

	process runRecalSNPApply {

  		tag "ALL"
  		// publishDir "${OUTDIR}/${params.tool}/Variants/Filtered"

  		input:
  		set file(recal_file),file(tranches),file(rscript),file(gvcf) from inputRecalSNPApply

  		output:
  		file vcf_snp   into outputRecalSNPApply

  		script:
 
  		vcf_snp = "genotypes.recal_SNP.vcf.gz"

  		"""
		gatk IndexFeatureFile -F $recal_file
	 	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
			-R $REF \
			-V $gvcf \
		        --recal-file $recal_file \
                	--tranches-file $tranches \
			-mode SNP \
			--ts-filter-level 99.0 \
			-O $vcf_snp	
  		"""
	}

	process runRecalIndelApply {

  		tag "ALL"
	  	// publishDir "${OUTDIR}/${params.tool}/Variants/Recal"

  		input:
	  	set file(recal_file),file(tranches),file(rscript),file(gvcf) from inputRecalIndelApply

	  	output:
	  	file vcf_indel into outputRecalIndelApply

	  	script:

  		vcf_indel = "genotypes.recal_Indel.vcf.gz"

  		"""
                gatk IndexFeatureFile -F $recal_file
        	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
                	-R $REF \
	                -V $gvcf \
        	        --recal-file $recal_file \
                	--tranches-file $tranches \
	                -mode INDEL \
        	        --ts-filter-level 99.0 \
                	-O $vcf_indel
	  	"""
	}

	process runVariantFiltrationIndel {

  		tag "ALL"
  		// publishDir "${OUTDIR}/${params.tool}/Variants/Filtered"

	  	input:
  		file(gvcf) from outputRecalIndelApply

	  	output:
	  	file(filtered_gvcf) into outputVariantFiltrationIndel

	  	script:

	  	filtered_gvcf = "genotypes.recal_Indel.filtered.vcf.gz"

  		"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration \
                -R $REF \
                -V $gvcf \
		--filter-expression "QD < 2.0" \
		--filter-name "QDFilter" \
                -O $filtered_gvcf
  		"""
	}

	inputCombineVariants = outputVariantFiltrationIndel.mix(outputRecalSNPApply)

	process runCombineVariants {

  		tag "ALL"
	  	// publishDir "${OUTDIR}/${params.tool}/Variants/Final", mode: 'copy'

  		input: 
	     	set file(indel),file(snp) from inputCombineVariants.collect()

  		output:
		set file(merged_file),file(merged_file_index) into inputLeftNormalize

		script:
		merged_file = "merged_callset.vcf.gz"
		merged_file_index = merged_file + ".tbi"

		"""

			gatk SortVcf -I $indel -O indels.sorted.vcf.gz
			gatk SortVcf -I $snp -O snps.sorted.vcf.gz
			gatk MergeVcfs \
			-I=indels.sorted.vcf.gz \
			-I=snps.sorted.vcf.gz \
			-O=merged.vcf.gz \
			-R=$REF \

			gatk SelectVariants \
			-R $REF \
			-V merged.vcf.gz \
			-O $merged_file \
			--remove-unused-alternates true \
			--exclude-non-variants true
  		"""
	}

// ++++++++++++++++++
// GATK3 workflow -  for legacy purposes only
// ++++++++++++++++++

} else if (params.tool == "gatk3") {

   process runBaseRecalibratorGATK3 {

    	tag "${indivID}|${sampleID}"
    	// publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Processing/BaseRecalibrator/", mode: 'copy'
	    
    	input:
    	set indivID, sampleID, dedup_bam, dedup_bai from MarkDuplicatesOutput
    
    	output:
    	set indivID, sampleID, dedup_bam, file(recal_table) into runBaseRecalibratorOutput
    
    	script:
    	recal_table = sampleID + "_recal_table.txt" 
       
    	"""
		java -XX:ParallelGCThreads=2 -Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=tmp/ -jar ${GATK} \
			-T BaseRecalibrator \
			-R ${REF} \
			-I ${dedup_bam} \
			-knownSites ${GOLD1} \
			-knownSites ${DBSNP} \
	                -knownSites ${G1K} \
			-o ${recal_table} \
			-nct ${task.cpus}
		"""
    }

    process runPrintReadsGATK3 {

    	tag "${indivID}|${sampleID}"
    	publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/", mode: 'copy'

    	scratch use_scratch
	    
    	input:
    	set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput 

    	output:
    	set indivID, sampleID, file(outfile_bam), file(outfile_bai) into runPrintReadsOutput_for_DepthOfCoverage, runPrintReadsOutput_for_Multiple_Metrics, inputHCSample
    	set indivID, sampleID, realign_bam, recal_table into runPrintReadsOutput_for_PostRecal
            
    	script:
	outfile_bam = sampleID + ".clean.bam"
    	outfile_bai = sampleID + ".clean.bai"
           
    	"""

		java -XX:ParallelGCThreads=1 -Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T PrintReads \
		-R ${REF} \
		-I ${realign_bam} \
		-BQSR ${recal_table} \
		-o ${outfile_bam} \
		-nct ${task.cpus}
    	"""
    }

    process runBaseRecalibratorPostRecalGATK3 {

    	tag "${indivID}|${sampleID}"
    	publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Processing/BaseRecalibratorPostRecal/", mode: 'copy'
	    
    	input:
    	set indivID, sampleID, realign_bam, recal_table from runPrintReadsOutput_for_PostRecal
    
    	output:
    	set indivID, sampleID, recal_table, file(post_recal_table) into runBaseRecalibratorPostRecalOutput_Analyze
        
    	script:
    	post_recal_table = sampleID + "_post_recal_table.txt" 
   
    	"""
    
	java -XX:ParallelGCThreads=1 -Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${DBSNP} \
		-BQSR ${recal_table} \
		-o ${post_recal_table}
	"""
    }
    process runAnalyzeCovariatesGATK3 {
	tag "${indivID}|${sampleID}"
    	publishDir "${OUTDIR}/${params.tool}/${indivID}/${sampleID}/Processing/AnalyzeCovariates/", mode: 'copy'
	    
    	input:
    	set indivID, sampleID, recal_table, post_recal_table from runBaseRecalibratorPostRecalOutput_Analyze

	output:
	set indivID, sampleID, recal_plots into runAnalyzeCovariatesOutput
	    
    	script:
    	recal_plots = sampleID + "_recal_plots.pdf" 

    	"""
    
	java -XX:ParallelGCThreads=1 -Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T AnalyzeCovariates \
		-R ${REF} \
		-before ${recal_table} \
		-after ${post_recal_table} \
		-plots ${recal_plots}
    	"""
    }

    process runHCSampleGATK3 {

  	tag "${SampleID}"
  	publishDir "${OUTDIR}/${params.tool}/${IndivID}/${SampleID}/HaplotypeCaller/" , mode: 'copy'

  	input: 
  	set IndivID,SampleID,file(bam),file(bai) from inputHCSample

  	output:
  	file(vcf) into outputHCSample

  	script:
  
  	vcf = SampleID + ".raw_variants.g.vcf"

  	"""
		java -jar -Xmx${task.memory.toGiga()}G $GATK \
		-T HaplotypeCaller \
		-R $REF \
		-I $bam \
		-L $TARGETS \
		-L chrM \
		--genotyping_mode DISCOVERY \
		--emitRefConfidence GVCF \
    		-o $vcf \
		-nct ${task.cpus}
  	"""
    }

    inputHCJoined = outputHCSample.collect()

    process runJoinedGenotypingGATK3 {
  
  	tag "ALL - using 17 IKMB reference exomes for calibration"
  	publishDir "${OUTDIR}/gatk3/JoinedGenotypes"
  
  	input:
  	file(vcf_list) from inputHCJoined
  
  	output:
  	file(gvcf) into (inputRecalSNP , inputRecalIndel)
  
  	script:
  
  	gvcf = "genotypes.gvcf"
  
  	"""
	 java -jar -Xmx${task.memory.toGiga()}G $GATK \
                -T GenotypeGVCFs \
                -R $REF \
                --variant ${vcf_list.join(" --variant ")} \
		--variant $calibration_exomes \
		--dbsnp $DBSNP \
                -o $gvcf \
		-nt ${task.cpus} \
		--useNewAFCalculator
  	"""
    }

    process runRecalibrationModeSNPGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Recal"

  	input:
  	file(gvcf) from inputRecalSNP

  	output:
  	set file(recal_file),file(tranches),file(rscript),file(gvcf) into inputRecalSNPApply

  	script:
  	recal_file = "genotypes.recal_SNP.recal"
  	tranches = "genotypes.recal_SNP.tranches"
  	rscript = "genotypes.recal_SNP.R"

  	"""
		java -jar -Xmx${task.memory.toGiga()}G $GATK -T SelectVariants \
			-R $REF \
			-V $gvcf \
			-o snps.vcf \
			-selectType SNP

		java -jar -Xmx${task.memory.toGiga()}G $GATK -T VariantRecalibrator \
		-R $REF \
		-input snps.vcf \
                --recal_file $recal_file \
                --tranches_file $tranches \
                --rscript_file $rscript \
		-an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP -an QD \
                --mode SNP \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
		-resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1K \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
		-nt ${task.cpus}
  	"""
     }

    process runRecalibrationModeIndelGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Recal"

  	input:
  	file(gvcf) from inputRecalIndel

  	output:
  	set file(recal_file),file(tranches),file(rscript),file(gvcf) into inputRecalIndelApply

  	script:

  	recal_file = "genotypes.recal_Indel.recal"
  	tranches = "genotypes.recal_Indel.tranches"
  	rscript = "genotypes.recal_Indel.R"

  	"""

                java -jar -Xmx${task.memory.toGiga()}G $GATK -T SelectVariants \
                        -R $REF \
                        -V $gvcf \
                        -o snps.vcf \
                        -selectType INDEL \
			-selectType MIXED \
			-selectType MNP \
			-selectType SYMBOLIC \
			-selectType NO_VARIATION

        	java -jar -Xmx${task.memory.toGiga()}G $GATK -T VariantRecalibrator \
                -R $REF \
                -input $gvcf \
                --recal_file $recal_file \
                --tranches_file $tranches \
                --rscript_file $rscript \
                -an MQ -an MQRankSum -an SOR -an ReadPosRankSum -an FS  \
                --mode INDEL \
                -resource:mills,known=false,training=true,truth=true,prior=15.0 $GOLD1 \
                -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
                -nt ${task.cpus}
  	"""
    }

    process runRecalSNPApplyGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Filtered"

  	input:
  	set file(recal_file),file(tranches),file(rscript),file(gvcf) from inputRecalSNPApply

  	output:
  	file vcf_snp   into outputRecalSNPApply

  	script:
 
  	vcf_snp = "genotypes.recal_SNP.vcf.gz"

  	"""
	 java -jar -Xmx${task.memory.toGiga()}G $GATK -T ApplyRecalibration \
		-R $REF \
		-input $gvcf \
	        -recalFile $recal_file \
                -tranchesFile $tranches \
		--mode SNP \
		--ts_filter_level 99.0 \
		-o $vcf_snp	
  	"""
    }

    process runRecalIndelApplyGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Recal"

  	input:
  	set file(recal_file),file(tranches),file(rscript),file(gvcf) from inputRecalIndelApply

  	output:
  	set file(vcf_indel),file(vcf_index) into outputRecalIndelApply

  	script:

  	vcf_indel = "genotypes.recal_Indel.vcf.gz"
	vcf_index = vcf_indel + ".tbi"

  	"""
         java -jar -Xmx${task.memory.toGiga()}G $GATK -T ApplyRecalibration \
                -R $REF \
                -input $gvcf \
                -recalFile $recal_file \
                -tranchesFile $tranches \
                --mode Indel \
                --ts_filter_level 99.0 \
                -o $vcf_indel
  	"""
    }

    process runVariantFiltrationIndelGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Filtered"

  	input:
  	set file(gvcf),file(idx) from outputRecalIndelApply

  	output:
  	file(filtered_gvcf) into outputVariantFiltrationIndel

  	script:

  	filtered_gvcf = "genotypes.recal_Indel.filtered.vcf.gz"

  	"""
		java -jar -Xmx${task.memory.toGiga()}G $GATK -T VariantFiltration \
                -R $REF \
                -V $gvcf \
		--filterName GATKStandardQD \
                --filterExpression "QD < 2.0" \
                --filterName GATKStandardReadPosRankSum \
                --filterExpression "ReadPosRankSum < -20.0" \
                --filterName GATKStandardFS \
                --filterExpression "FS > 200.0" \
                -o $filtered_gvcf
  	"""
    }

    inputCombineVariants = outputVariantFiltrationIndel.mix(outputRecalSNPApply)

    process runCombineVariantsGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Final", mode: 'copy'

  	input: 
     	set file(indel),file(snp) from inputCombineVariants.collect()

  	output:
    	set file(merged_file),file(merged_file_index) into outputCombineVariants

  	script:
    	merged_file = "merged_callset.vcf.gz"
	merged_file_index = merged_file + ".tbi"

  	"""
		tabix $indel
		tabix $snp
		
		java -jar -Xmx${task.memory.toGiga()}G $GATK -T CombineVariants \
		-R $REF \
		--variant $indel --variant $snp \
		-o $merged_file \
		--genotypemergeoption UNSORTED
  	"""
    }

    process runRemoveCalibrationExomesGATK3 {

  	tag "ALL"
  	// publishDir "${OUTDIR}/${params.tool}/Final", mode: 'copy'

  	input:
  	set file(merged_vcf),file(merged_vcf_index) from outputCombineVariants

  	output:
  	set file(filtered_vcf),file(filtered_vcf_index) into inputLeftNormalize

  	script:
  	filtered_vcf = "merged_callset.calibration_removed.vcf.gz"
	filtered_vcf_index = filtered_vcf + ".tbi"

  	"""
		java -jar -Xmx${task.memory.toGiga()}G $GATK -T SelectVariants \
        	-R $REF \
		-V $merged_vcf \
		-L chrM \
		-L $TARGETS \
		--exclude_sample_file $calibration_samples_list \
        	-o $filtered_vcf \
		-env \
		-trimAlternates
  	"""
    }

} 

// *********************
// Compute statistics for fastQ files, libraries and samples
// *********************

process runCollectMultipleMetrics {
	tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'
 
	scratch use_scratch
	    
	input:
	set indivID, sampleID, bam, bai from BamForMultipleMetrics

	output:
	file("${prefix}*") into CollectMultipleMetricsOutput mode flatten

	script:       
	prefix = sampleID + "."

	"""
		java -XX:ParallelGCThreads=1 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
       	        PROGRAM=CollectSequencingArtifactMetrics \
                PROGRAM=CollectQualityYieldMetrics \
	        PROGRAM=CollectGcBiasMetrics \
		PROGRAM=CollectBaseDistributionByCycle \
		INPUT=${bam} \
		REFERENCE_SEQUENCE=${REF} \
		DB_SNP=${DBSNP} \
		INTERVALS=${BAITS} \
		ASSUME_SORTED=true \
		QUIET=true \
		OUTPUT=${prefix} \
		TMP_DIR=tmp
	"""
}	

process runHybridCaptureMetrics {

    tag "${indivID}|${sampleID}"
    publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_HC_Metrics

    output:
    file(outfile) into HybridCaptureMetricsOutput mode flatten

    script:
    outfile = sampleID + ".hybrid_selection_metrics.txt"

    """
        java -XX:ParallelGCThreads=1 -Xmx10g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectHsMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                TARGET_INTERVALS=${TARGETS} \
                BAIT_INTERVALS=${BAITS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

process runOxoGMetrics {

    tag "${indivID}|${sampleID}"
    publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_OxoG_Metrics

    output:
    file(outfile) into runOxoGMetricsOutput mode flatten

    script:
    outfile = sampleID + ".OxoG_metrics.txt"

    """

        java -XX:ParallelGCThreads=1 -Xmx10g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectOxoGMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                DB_SNP=${DBSNP} \
                INTERVALS=${BAITS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

process runDepthOfCoverage {

       tag "${indivID}|${sampleID}"
       publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/DepthOfCoverage", mode: 'copy'

       scratch use_scratch

       input:
       set indivID, sampleID, file(bam), file(bai) from BamForDepthOfCoverage

       output:
       file("${prefix}*") into DepthOfCoverageOutput

       script:
       prefix = sampleID + "."

       """

       java -XX:ParallelGCThreads=1 -Djava.io.tmpdir=tmp/ -Xmx10g -jar ${GATK} \
             -R ${REF} \
             -T DepthOfCoverage \
             -I ${bam} \
             --omitDepthOutputAtEachBase \
             -ct 10 -ct 20 -ct 50 -ct 100 \
             -o ${sampleID}

       """
}

process runFastQC {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/Common/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastQC/", mode: 'copy'

    scratch use_scratch

    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastQC

    output:
    set file("*.zip"), file("*.html") into FastQCOutput
    	
    script:
    """
    fastqc -t 1 -o . ${fastqR1} ${fastqR2}
    """
}



// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

process runMultiQCFastq {

    tag "Generating fastq level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Fastq", mode: 'copy'
	    
    input:
    file('*') from FastQCOutput.flatten().toList()
    
    output:
    file("fastq_multiqc*") into runMultiQCFastqOutput
    	
    script:

    """
   cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n fastq_multiqc *.zip *.html
    """
}

process runMultiQCLibrary {

    tag "Generating library level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Library", mode: 'copy'
	    
    input:
    file('*') from DuplicatesOutput_QC.flatten().toList()

    output:
    file("library_multiqc*") into runMultiQCLibraryOutput
    	
    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n library_multiqc *.txt
    """
}

process runMultiQCSample {

    tag "Generating sample level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Sample", mode: 'copy'
	    
    input:
    file('*') from CollectMultipleMetricsOutput.flatten().toList()
    file('*') from HybridCaptureMetricsOutput.flatten().toList()
    file('*') from runOxoGMetricsOutput.flatten().toList()
        
    output:
    file("sample_multiqc*") into runMultiQCSampleOutput
    	
    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n sample_multiqc *
    """
}


// *************************
// Variant effect prediction
// *************************

// Left-normalize the variants
process runLeftNormalize {

   tag "ALL"
   publishDir "${OUTDIR}/${params.tool}/Final", mode: 'copy'

   input:
   set file(vcf_file),file(index) from inputLeftNormalize
 
   output:
   set file(vcf_normalized),file(vcf_normalized_index) into ( inputVep, inputAnnovar)

   script:

   vcf_normalized = "variants.merged.normalized.vcf.gz"
   vcf_normalized_index = vcf_normalized + ".tbi"

   """
	bcftools norm -f $REF $vcf_file | bgzip > $vcf_normalized
	tabix $vcf_normalized
   """

}

process runVep {

 tag "ALL"
 publishDir "${OUTDIR}/${params.tool}/Annotation/VEP", mode: 'copy'
 
 input:
   set file(vcf_file),file(vcf_index) from inputVep

 output:
   file('annotation.vep.vcf') into outputVep

 script:

   """
     vep --offline --cache --dir $VEP_CACHE --fork ${task.cpus} \
 	--assembly GRCh37 -i $vcf_file -o annotation.vep.vcf --allele_number --canonical \
	--force_overwrite --vcf --no_progress \
	--pubmed \
	--plugin ExAC,$EXAC \
	--plugin CADD,$CADD \
	--plugin LoFtool --plugin LoF \
	--fasta /ifs/data/nfs_share/ikmb_repository/references/genomes/homo_sapiens/EnsEMBL/GRCh37/genome.fa
   """

}

process runAnnovar {

 tag "ALL"
 publishDir "${OUTDIR}/${params.tool}/Annotation/Annovar", mode: 'copy'

 input:
   set file(vcf_file),file(vcf_index) from inputAnnovar

 output:
   file(annovar_result) into outputAnnovar

 script:
  annovar_target = vcf_file + ".annovar"
  annovar_result = vcf_file + ".annovar.hg19_multianno.vcf"

   """
      table_annovar.pl -v \
	--protocol ensGene,knownGene,refGene,dbnsfp33a,intervar_20170202,esp6500siv2_all,gnomad_exome,clinvar_20170905,cadd13gt10 \
	--operation g,g,g,f,f,f,f,f,f \
	--outfile $annovar_target \
	--buildver hg19 \
	--remove \
        --thread ${task.cpus} \
	--otherinfo \
	--vcfinput \
	${vcf_file} $ANNOVAR_DB
   """

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  if (params.email) {

            def subject = 'Diagnostic exome analysis finished.'
            def recipient = params.email

            ['mail', '-s', subject, recipient].execute() << """

            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
            """
  }

}


//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

