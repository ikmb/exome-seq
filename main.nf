#!/usr/bin/env nextflow

/**
===============================
IKMB Diagnostic Exome Pipeline
===============================

This Pipeline performs one of two workflows to generate variant calls and effect predictions
using either the GATK processing chain or Freebayes.

### Homepage / git
http://git.ikmb.uni-kiel.de/bfx-core/NF-diagnostics-exome
### Implementation
Implemented in Q1 2018

This pipeline is based on the updated GATK best-practices (where applicable).
 - trimming (Trimgalore)
 - Alignment (BWA)
 - Duplicate marking (GATK)
 - recalibration 
 - variant calling
 - variant recalibration and filtering
 - variant effect prediction

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Pipeline version
VERSION = "1.0"
params.version = VERSION

// Help message
helpMessage = """
===============================================================================
IKMB Diagnostic Exome pipeline | version ${VERSION}
===============================================================================
Usage: nextflow -c /path/to/git/nextflow.config run /path/to/git/main.nf --assembly hg19_clinical --kit Nextera --samples Samples.csv
This example will perform an exome analysis against the hg19 (with decoys) assembly, assuming that exome reads were generated with
the Nextera kit and using the GATK4 best-practice workflow. 
Required parameters:
--samples                      A sample list in CSV format (see website for formatting hints)
--assembly                     Name of the reference assembly to use
--kit			       Name of the exome kit (available options: xGen, xGen_custom, Nextera)
Optional parameters:
--effect_prediction            Whether to run effect prediction on the final variant set (default: false)
--vqsr 			       Whether to also run variant score recalibration (only works >= 30 samples) (default: false)
--run_name 		       A descriptive name for this pipeline run
--cram			       Whether to output the alignments in CRAM format (default: bam)
--fasta			       A reference genome in FASTA format (set automatically if using --assembly)
--dbsnp			       dbSNP data in VCF format (set automatically if using --assembly)
--g1k			       A SNP reference (usually 1000genomes, set automatically if using --assembly)
--mills_indels		       An INDEL reference (usually MILLS/1000genomes, set automatically if using --assembly)
--omni			       An SNP reference (usually OMNI, set automatically if using --assembly)
--hapmap		       A SNP reference (usually HAPMAP, set automatically if using --assembly)
--targets		       A interval_list target file (set automatically if using the --kit option)
--baits			       A interval_list bait file (set automatically if using the --kit option)
--interval_padding	       For GATK, include this number of nt upstream and downstream around the exome targets (default: 50)
Output:
--outdir                       Local directory to which all output is written (default: output)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// #############
// INPUT OPTIONS
// #############

// Sample input file
inputFile = file(params.samples)

// Giving this pipeline run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

// This will eventually enable switching between multiple assembly versions
// Currently, only hg19 has all the required reference files available
params.assembly = "hg19"

REF = params.fasta ?: file(params.genomes[ params.assembly ].fasta)
DBSNP = params.dbsnp ?: file(params.genomes[ params.assembly ].dbsnp )
G1K = params.g1k ?: file(params.genomes[ params.assembly ].g1k )
MILLS = params.mills_indels ?: file(params.genomes[ params.assembly ].mills )
OMNI = params.omni ?: file(params.genomes[ params.assembly ].omni )
HAPMAP = params.hapmap ?: file(params.genomes[ params.assembly ].hapmap )
VEP_CACHE = params.vep_cache
MITOCHONDRION = params.mitochondrion ?: params.genomes[ params.assembly ].mitochondrion

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

SNP_RULES = params.snp_filter_rules
INDEL_RULES = params.indel_filter_rules

// Annotations to use for variant recalibration
snp_recalibration_values = params.snp_recalibration_values
indel_recalbration_values = params.indel_recalbration_values

params.effect_prediction = true
params.hard_filter = false

// Whether to produce BAM output instead of CRAM
params.cram = false
align_suffix = (params.cram == false) ? "bam" : "cram"

// Location of applications used
OUTDIR = file(params.outdir)

// Available exome kits

if (TARGETS == false || BAITS == false ) {
   exit 1, "Information on enrichment kit incomplete or missing (please see the documentation for details!"
}

// Whether to send a notification upon workflow completion
params.email = false

if(params.email == false) {
	exit 1, "You must provide an Email address to which pipeline updates are send!"
}

// Whether to use a local scratch disc
use_scratch = params.scratch

// Make sure the Nextflow version is current enough
try {
    if( ! nextflow.version.matches(">= $params.nextflow_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nextflow_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please use a more recent version of Nextflow!\n" +
              "============================================================"
}

logParams(params, "pipeline_parameters.txt")

// Header log info
log.info "========================================="
log.info "IKMB Diagnostic Exome pipeline v${VERSION}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${params.run_name}"
log.info "========================================="

// Read sample file 
Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  readPairsFastp }

process runFastp {

	tag "${indivID}|${sampleID}"

	input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, date, fastqR1, fastqR2 from readPairsFastp

	output:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, date, center, file(left),file(right) into inputBwa
	set file(html),file(json) into fastp_results

	script:
	left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
	json = file(fastqR1).getBaseName() + ".fastp.json"
	html = file(fastqR1).getBaseName() + ".fastp.html"

	"""
		fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runBWA {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    // publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/", mode: 'copy'

    // scratch use_scratch
	
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa
    
    output:
    set indivID, sampleID, file(outfile) into runBWAOutput
    
    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"	

    """
	bwa mem -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t ${task.cpus} ${REF} $left $right | samtools sort -O bam -m 2G -@ 4 - > $outfile
    """	
}

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])

process mergeBamFiles_bySample {

        tag "${indivID}|${sampleID}"
	
	input:
        set indivID, sampleID, file(aligned_bam_list) from runBWAOutput_grouped_by_sample

	output:
	set indivID,sampleID,file(merged_bam),file(merged_bam_index) into mergedBamFile_by_Sample

	script:
	merged_bam = sampleID + ".merged.bam"
	merged_bam_index = merged_bam + ".bai"

	"""

	    	gatk MergeSamFiles \
                    -I ${aligned_bam_list.join(' -I ')} \
                    -O merged.bam \
		    --USE_THREADING true \
                    --SORT_ORDER coordinate

		gatk SetNmMdAndUqTags \
			-I merged.bam \
			-O $merged_bam \
			-R $REF \
			--IS_BISULFITE_SEQUENCE false

		samtools index $merged_bam

	"""
}

process runMarkDuplicates {

	tag "${indivID}|${sampleID}"
        publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates", mode: 'copy'

        scratch use_scratch

        input:
        set indivID, sampleID, file(merged_bam),file(merged_bam_index) from mergedBamFile_by_Sample

        output:
        set indivID, sampleID, file(outfile_bam),file(outfile_bai) into MarkDuplicatesOutput, BamForMultipleMetrics, runPrintReadsOutput_for_OxoG_Metrics, runPrintReadsOutput_for_HC_Metrics, BamForDepthOfCoverage
	file(outfile_md5) into MarkDuplicatesMD5
	file(outfile_metrics) into DuplicatesOutput_QC

        script:
        outfile_bam = sampleID + ".dedup.bam"
        outfile_bai = sampleID + ".dedup.bai"
	outfile_md5 = sampleID + ".dedup.bam.md5"

        outfile_metrics = sampleID + "_duplicate_metrics.txt"

	"""
        	gatk --java-options "-Xmx${task.memory.toGiga()-1}G" MarkDuplicates \
                	-I ${merged_bam} \
	                -O ${outfile_bam} \
        	        -M ${outfile_metrics} \
                        --CREATE_INDEX true \
			--ASSUME_SORT_ORDER=coordinate \
			--MAX_RECORDS_IN_RAM 100000 \
			--CREATE_MD5_FILE true \
                        --TMP_DIR tmp \
			-R ${REF}
	"""

}

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
	// publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibrator/", mode: 'copy'
	    
	input:
	set indivID, sampleID, dedup_bam, dedup_bai from MarkDuplicatesOutput
    
	output:
	set indivID, sampleID, dedup_bam, file(recal_table) into runBaseRecalibratorOutput

	script:
	recal_table = sampleID + "_recal_table.txt" 

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" BaseRecalibrator \
		--reference ${REF} \
		-L $TARGETS \
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
		--use-original-qualities \
		--input ${dedup_bam} \
		--known-sites ${MILLS} \
		--known-sites ${DBSNP} \
		--output ${recal_table}
	"""
}

process runApplyBQSR {

	tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/", mode: 'copy'

	scratch use_scratch
	    
	input:
	set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput 

	output:
	set indivID, sampleID, file(outfile_bam), file("*.bai") into runPrintReadsOutput_for_Multiple_Metrics,inputHCSample,inputCollectReadCounts
	set indivID, sampleID, realign_bam, recal_table into runPrintReadsOutput_for_PostRecal
	set indivID, outfile_md5 into BamMD5
            
	script:

	outfile_bam = sampleID + ".clean.${align_suffix}"
	outfile_bai = sampleID + ".clean.${align_suffix}.bai"
	outfile_md5 = sampleID + ".clean.${align_suffix}.md5"
           
    	"""
        	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyBQSR \
                --reference ${REF} \
		--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
                --input ${realign_bam} \
		--use-original-qualities \
		-OBI true \
		-L $TARGETS \
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
                -bqsr ${recal_table} \
                --output ${outfile_bam} \
                -OBM true \
    	"""
}    

// Call variants on a per-sample basis

process runHCSample {

	tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Variants/HaplotypeCaller" , mode: 'copy'

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
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
		--emit-ref-confidence GVCF \
		-OVI true \
    		--output $vcf \
		--native-pair-hmm-threads ${task.cpus} &> log.txt \
  	"""
}

// Import individual vcf files into a GenomicsDB database on a per chromosome basis
// From here on all samples are in the same file
process runGenomicsDBImport  {

	tag "ALL"
        publishDir "${OUTDIR}/GATK/Variants/JointGenotypes/", mode: 'copy'

	scratch use_scratch 

	input:
        file(vcf_list) from outputHCSample.collect()
	file(index_list) from outputHCSampleIndex.collect()

	output:
        set file(merged_vcf),file(merged_vcf_index) into inputGenotypeGVCFs

	script:
	merged_vcf = "merged.g.vcf.gz"
	merged_vcf_index = merged_vcf + ".tbi"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" CombineGVCFs  \
		--variant ${vcf_list.join(" --variant ")} \
		--reference $REF \
		--intervals $TARGETS \
		--intervals $MITOCHONDRION \
		-ip ${params.interval_padding} \
		--OVI true \
		--output $merged_vcf \
	"""

}

// Perform genotyping on a per chromosome basis

process runGenotypeGVCFs {
  
	tag "ALL"
	publishDir "${OUTDIR}/GATK/Variants/JointGenotypes", mode: 'copy'
  
	input:
	set file(merged_vcf), file(merged_vcf_index) from inputGenotypeGVCFs
  
	output:
	set file(gvcf), file(gvcf_index) into (inputHardFilterSNP, inputRecalSNP, inputHardFilterIndel, inputRecalIndel )
  
	script:
  
	gvcf = "genotypes.vcf.gz"
	gvcf_index = gvcf + ".tbi"
  
	"""
 	gatk --java-options "-Xmx${task.memory.toGiga()}G" GenotypeGVCFs \
		--reference $REF \
		--dbsnp $DBSNP \
		-L $TARGETS \
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
		-new-qual \
		--only-output-calls-starting-in-intervals \
		-V $merged_vcf \
              	--output $gvcf \
                -G StandardAnnotation \
		-OVI true
	"""
}

////////////////////////
// Hard filtering
////////////////////////

process runHardFilterSNP {
		
	tag "ALL"
	publishDir "${OUTDIR}/GATK/Variants/HardFilter/", mode: 'copy'

	input:
	set file(vcf),file(vcf_index) from inputHardFilterSNP

	output:
	set file(vcf_filtered),file(vcf_filtered_index) into outputHardFilterSNP

	script:
	vcf_filtered = "genotypes.merged.snps.filtered.vcf.gz"
	vcf_filtered_index = vcf_filtered + ".tbi"

	"""
		gatk SelectVariants \
			-R $REF \
			-V $vcf \
			--select-type-to-include SNP \
			--output genotypes.merged.snps.vcf.gz \
			-OVI true
				
		gatk VariantFiltration \
			-R $REF \
			-V genotypes.merged.snps.vcf.gz \
			--output $vcf_filtered \
			--filter-expression "${SNP_RULES}" \
			--filter-name "hard_snp_filter" \
			-OVI true
	"""

}

process runHardFilterIndel {

	tag "ALL"
        publishDir "${OUTDIR}/GATK/Variants/HardFilter", mode: 'copy'
        
        input:
        set file(vcf),file(vcf_index) from inputHardFilterIndel

        output:
        set file(vcf_filtered),file(vcf_filtered_index) into outputHardFilterIndel

        script:
        vcf_filtered = "genotypes.merged.indels.filtered.vcf.gz"
        vcf_filtered_index = vcf_filtered + ".tbi"

        """
 	       gatk SelectVariants \
               -R $REF \
               -V $vcf \
               --select-type-to-include INDEL \
               --output genotypes.merged.indels.vcf.gz \
               -OVI true

              gatk VariantFiltration \
              -R $REF \
              -V genotypes.merged.indels.vcf.gz \
              --output $vcf_filtered \
	      --filter-expression "${INDEL_RULES}" \
              --filter-name "hard_indel_filter" \
	      -OVI true
        """
}

process runCombineHardVariants {

	tag "ALL"
        publishDir "${OUTDIR}/GATK/Variants/HardFilter/Preprocess", mode: 'copy'

        input:
        set file(indel),file(indel_index) from outputHardFilterIndel
	set file(snp),file(snp_index) from outputHardFilterSNP

        output:
        set file(merged_file),file(merged_file_index) into inputfilterPassVariants

        script:
        merged_file = "${run_name}.merged_callset.hard.vcf.gz"
        merged_file_index = merged_file + ".tbi"

        """
		gatk --java-options "-Xmx${task.memory.toGiga()}G"  MergeVcfs \
		-I $indel \
		-I $snp \
		-O $merged_file \
		
		gatk IndexFeatureFile -F $merged_file

        """
}

process runFilterPassVariants {

	tag "ALL"
        publishDir "${OUTDIR}/GATK/Variants/HardFilter/Preprocess", mode: 'copy'

	input:
	set file(vcf),file(index) from inputfilterPassVariants

	output:
	set file(vcf_pass),file(vcf_pass_index) into inputSplitHardVariants

	script:
	vcf_pass = "${run_name}.merged_callset.hard.pass.vcf.gz"
	vcf_pass_index = vcf_pass + ".tbi"

	"""
		bcftools view -f "PASS" $vcf | bgzip > $vcf_pass
		tabix $vcf_pass
	"""

}

process runSplitHardVariantsBySample {

	tag "ALL|${params.assembly}"
        publishDir "${OUTDIR}/GATK/Variants/HardFilter/Final/BySample", mode: 'copy'

        input:
        set file(vcf_clean),file(vcf_clean_index) from inputSplitHardVariants

        output:
        file("*.vcf.gz") into HardVcfBySample

        script:

        """
                for sample in `bcftools query -l ${vcf_clean}`; do gatk SelectVariants -R $REF -V ${vcf_clean} --exclude-non-variants --remove-unused-alternates -sn \$sample -O \$sample'.vcf.gz' ; done;
        """

}


/////////////////////////
// Variant recalibration 
/////////////////////////

process runRecalibrationModeSNP {

	tag "ALL"
	publishDir "${OUTDIR}/GATK/Variants/Recal"
	
	input:
	set file(vcf),file(vcf_index) from inputRecalSNP

	output:
	set file(recal_file),file(tranches) into inputRecalSNPApply

	when:
	params.vqsr == true

	script:
	recal_file = "genotypes.merged.snps.recal"
	tranches = "genotypes.merged.snps.tranches"

	"""

		gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
		-R $REF \
		-V $vcf \
	       	-O $recal_file \
       	        --tranches-file $tranches \
		-an ${snp_recalibration_values.join(' -an ')} \
	        -mode SNP \
		-OVI true \
		--resource hapmap,known=false,training=true,truth=true,prior=15.0:$HAPMAP \
		--resource omni,known=false,training=true,truth=true,prior=12.0:$OMNI \
		--resource 1000G,known=false,training=true,truth=false,prior=10.0:$G1K \
		--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
		-tranche ${params.snp_recalibration_tranche_values.join(' -tranche ')} \
		--max-gaussians 4
	"""
}

process runRecalibrationModeIndel {
	
	tag "ALL"
	publishDir "${OUTDIR}/GATK/Variants/Recal"

  	input:
	set file(vcf),file(vcf_index) from inputRecalIndel

  	output:
	set file(recal_file),file(tranches),file(vcf),file(vcf_index) into inputRecalIndelApply

	when:
	params.vqsr == true

	script:
  	recal_file = "genotypes.merged.indel.recal"
	tranches = "genotypes.merged.indel.tranches"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
        	        -R $REF \
	        	-V $vcf \
	               	-O $recal_file \
        	        --tranches-file $tranches \
	                -an ${indel_recalibration_values.join(' -an ')} \
	                -mode INDEL \
			-OVI true \
	        	--resource mills,known=false,training=true,truth=true,prior=15.0:$MILLS \
	               	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
			-tranche ${params.indel_recalibration_tranche_values.join(' -tranche ')} \
			--max-gaussians 3
  	"""
}

process runRecalIndelApply {

	tag "ALL"
        publishDir "${OUTDIR}/GATK/Variants/Recal"

        input:
        set file(recal_file),file(tranches),file(gvcf),file(gvcf_index) from inputRecalIndelApply

        output:
        set file(vcf_indel),file(vcf_indel_index) into outputRecalIndelApply

        script:

        vcf_indel = "genotypes.recal_Indel.vcf.gz"
        vcf_indel_index = vcf_indel + ".tbi"

        """
        	gatk IndexFeatureFile -F $recal_file
                        gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
                        -R $REF \
                        -V $gvcf \
                        --recal-file $recal_file \
                        --tranches-file $tranches \
                        -mode INDEL \
                        --ts-filter-level ${params.indel_filter_level} \
                        -OVI true \
                         -O $vcf_indel
        """
}

process runRecalSNPApply {
	
	tag "ALL"
	publishDir "${OUTDIR}/GATK/Variants/Filtered"
	
	input:
	set file(vcf),file(index) from outputRecalIndelApply
	set file(recal_file),file(tranches) from inputRecalSNPApply

	output:
	set file(vcf_snp),file(vcf_snp_index) into outputRecalSNPApply

	script:
 
	vcf_snp = "genotypes.recal_Indel.recal_SNP.vcf.gz"
	vcf_snp_index = vcf_snp + ".tbi"

	"""
	gatk IndexFeatureFile -F $recal_file
	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
		-R $REF \
		-V $vcf \
	        --recal-file $recal_file \
       	       	--tranches-file $tranches \
		-mode SNP \
		--ts-filter-level ${params.snp_filter_level} \
		-O $vcf_snp \
		-OVI true	
	"""
}

process runVariantFiltrationIndel {

	tag "ALL"
	publishDir "${OUTDIR}/GATK/Variants/Filtered"

  	input:
	set file(vcf),file(vcf_index) from outputRecalIndelApply

  	output:
  	set file(filtered_gvcf),file(filtered_gvcf_index) into inputSelectVariants

  	script:

  	filtered_gvcf = "genotypes.filtered.final.vcf.gz"
	filtered_gvcf_index = filtered_gvcf + ".tbi"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration \
       	       -R $REF \
               	-V $vcf \
		--filter-expression "QD < 2.0" \
		--filter-name "QDFilter" \
               	--output $filtered_gvcf \
		-OVI true
  	"""
}

process runSelectVariants {

	tag "ALL|${params.assembly}"
	publishDir "${OUTDIR}/GATK/Variants/Final", mode: 'copy'

	input:
	set file(vcf),file(vcf_index) from inputSelectVariants

	output:
	set file(vcf_clean),file(vcf_clean_index) into inputVep, inputSplitSample

	script:
	vcf_clean = "${run_name}.variants.merged.filtered.controls_removed.vcf.gz"
	vcf_clean_index = vcf_clean + ".tbi"

	"""
		gatk SelectVariants \
		-V $vcf \
		-R $REF \
		-O $vcf_clean \
		--remove-unused-alternates true \
		-OVI true \
		--exclude-non-variants true \
	"""

}

process runSplitBySample {

        tag "ALL|${params.assembly}"
        publishDir "${OUTDIR}/GATK/Variants/Final/BySample", mode: 'copy'

	input:
	set file(vcf_clean),file(vcf_clean_index) from inputSplitSample

	output: 
	file("*.vcf.gz*") into VcfBySample

	script: 

	"""
		for sample in `bcftools query -l ${vcf_clean}`; do gatk SelectVariants -R $REF -V ${vcf_clean} --exclude-non-variants --remove-unused-alternates -sn \$sample -O \$sample'.vcf.gz' -OVI true ; done;
	"""

}

// *********************
// Compute statistics for fastQ files, libraries and samples
// *********************

process runCollectMultipleMetrics {
	tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'
 
	scratch use_scratch
	    
	input:
	set indivID, sampleID, bam, bai from BamForMultipleMetrics

	output:
	file("${prefix}*") into CollectMultipleMetricsOutput mode flatten

	script:       
	prefix = sampleID + "."

	"""
		picard -Xmx5g CollectMultipleMetrics \
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
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_HC_Metrics

    output:
    file(outfile) into HybridCaptureMetricsOutput mode flatten
    file(outfile_per_target) into HsMetricsPerTarget

    script:
    outfile = sampleID + ".hybrid_selection_metrics.txt"
    outfile_per_target = sampleID + ".hybrid_selection_per_target_metrics.txt"

    """
        picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
		PER_TARGET_COVERAGE=${outfile_per_target} \
                TARGET_INTERVALS=${TARGETS} \
                BAIT_INTERVALS=${BAITS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

process runOxoGMetrics {

    tag "${indivID}|${sampleID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_OxoG_Metrics

    output:
    file(outfile) into runOxoGMetricsOutput mode flatten

    script:
    outfile = sampleID + ".OxoG_metrics.txt"

    """

         picard -Xmx${task.memory.toGiga()}G CollectOxoGMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                DB_SNP=${DBSNP} \
                INTERVALS=${TARGETS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

process runMultiqcFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastq", mode: 'copy'
	    
    input:
    file('*') from fastp_results.flatten().toList()
    
    output:
    file("fastp_multiqc*") into runMultiQCFastqOutput
    	
    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n fastp_multiqc *.json *.html
    """
}

process runMultiqcLibrary {

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

process runMultiqcSample {

    tag "Generating sample level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Sample", mode: 'copy'
	    
    input:
    file('*') from CollectMultipleMetricsOutput.flatten().toList()
    file('*') from HybridCaptureMetricsOutput.flatten().toList()
    file('*') from runOxoGMetricsOutput.flatten().toList()
        
    output:
    file("sample_multiqc.html") into runMultiQCSampleOutput
    	
    script:

    def subject = 'Diagnostic exome analysis quality report'
    def recipient = params.email

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n sample_multiqc *

    """
}

// *************************
// Variant effect prediction
// *************************

process runVep {

 tag "ALL"
 publishDir "${OUTDIR}/Annotation/VEP", mode: 'copy'
 
input:
   file(vcf_file) from inputVep

 output:
   file(annotated_vcf) into outputVep

 when:
 	params.effect_prediction == true

 script:
   annotated_vcf = run_name + ".annotation.vep.vcf"

   """
      vep --offline --cache --dir $VEP_CACHE --fork ${task.cpus} \
 	--assembly GRCh37 -i $vcf_file -o $annotated_vcf --allele_number --canonical \
	--force_overwrite --vcf --no_progress \
	--merged \
	--pubmed \
	--plugin LoFtool --plugin LoF \
	--fasta ${params.vep_fasta}
   """

}

process runSplitVEPBySample {

        tag "ALL|${params.assembly}"
        publishDir "${OUTDIR}/Annotation/VEP/BySample", mode: 'copy'

        input:
        file(vcf) from outputVep

        output:
        file("*.vcf") into VepVcfBySample

        script:

        """
                for sample in `bcftools query -l ${vcf}`; do gatk SelectVariants -R $REF -V ${vcf} -sn \$sample -O \$sample'.vcf' ; done;
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

