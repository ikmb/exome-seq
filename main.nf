#!/usr/bin/env nextflow

/**
===============================
Exome Pipeline
===============================

This Pipeline performs one of two workflows to generate variant calls 
using Google DeepVariant

### Homepage / git
git@github.com:ikmb/exome-seq.git
### Implementation
Implemented in Q1 2019

This pipeline is based on the DeepVariant best-practices (where applicable).
 - trimming (FastP)
 - Alignment (BWA)
 - Duplicate marking (Samtools)
 - Variant calling per-sample (Deepvariant)
 - Variant calling across samples (GLNexus)

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Pipeline version

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
IKMB Diagnostic Exome pipeline | version ${params.version}
===============================================================================
Usage: nextflow run ikmb/exome-seq --assembly GRCh38 --kit xGen_v2 --samples Samples.csv
This example will perform an exome analysis against the ALT-free hg38 assembly, assuming that exome reads were generated with
the IDT xGen v2 kit and using DeepVariant with GLNexus.

Required parameters:
--samples                      A sample list in CSV format (see website for formatting hints)
--assembly                     Name of the reference assembly to use
--kit			       Name of the exome kit (available options: xGen, xGen_custom, xGen_v2, Nextera, Pan_cancer)
--email 		       Email address to send reports to (enclosed in '')
Optional parameters:
--joint_calling		       Perform joint calling of all samples 
--skip_multiqc		       Don't attached MultiQC report to the email. 
--panel 		       Gene panel to check coverage of (valid options: cardio_dilatative, cardio_hypertrophic, cardio_non_compaction, eoIBD_25kb, imm_eoIBD_full, breast_cancer)
--all_panels 		       Run all gene panels defined for this assembly (none if no panel is defined!)
--panel_intervals	       Run a custom gene panel in interval list format (must have a matching sequence dictionary!)
--run_name 		       A descriptive name for this pipeline run
--cram			       Whether to output the alignments in CRAM format (default: bam)
--interval_padding	       Include this number of nt upstream and downstream around the exome targets (default: 10)
Expert options (usually not necessary to change!):
--fasta                        A reference genome in FASTA format (set automatically if using --assembly)
--dict                         A sequence dictionary matching --fasta (set automatically if using --assembly)
--dbsnp                        dbSNP data in VCF format (set automatically if using --assembly)
--targets                      A interval_list target file (set automatically if using the --kit option)
--baits                        A interval_list bait file (set automatically if using the --kit option)
--bed                          A list of calling intervals to be used by Deepvariant (default: exome kit targets will be converted to bed)
--max_length                   Cut reads down to this length (optional, default 0 = no trimming)
--min_mapq		       Minimum mapping quality to consider for general coverage analysis (default = 20)
--kill                         A list of known bad exons in the genome build that are ignored for panel coverage statistics (see documentation for details)
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

// #############
// INPUT OPTIONS
// #############

// Sample input file
inputFile = file(params.samples)

// Giving this pipeline run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

if (params.run_name == false) {
	log.info "No run name was specified, using ${run_name} instead"
}

// This will eventually enable switching between multiple assembly versions
// Currently, only hg19 has all the required reference files available
params.assembly = "GRCh38"

FASTA = file(params.genomes[ params.assembly ].fasta)
FAI_F = file(params.genomes[ params.assembly ].fai)
FASTAGZ_F = file(params.genomes[ params.assembly ].fastagz)
GZFAI_F = file(params.genomes[ params.assembly ].gzfai)
GZI_F = file(params.genomes[ params.assembly ].gzi)
DICT = file(params.genomes[ params.assembly ].dict)
DBSNP = file(params.genomes[ params.assembly ].dbsnp )

MITOCHONDRION = params.mitochondrion ?: params.genomes[ params.assembly ].mitochondrion

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

targets_to_bed = Channel.fromPath(TARGETS)

Channel.from(file(TARGETS))
	.into { TargetsToHS; TargetsToMetrics; TargetsToOxo }

Channel.from(file(BAITS))
	.into { BaitsToHS; BaitsToMetrics }

if (params.kill) {
	KILL = params.kill
} else if (params.kit && params.genomes[params.assembly].kits[params.kit].kill) {
	KILL = params.genomes[params.assembly].kits[params.kit].kill
} else {
	KILL = false
}

/*
PANEL COVERAGE - pick the correct panel for reporting
*/

if (params.panel && params.all_panels || params.panel && params.panel_intervals || params.all_panels && params.panel_intervals) {
	log.info "The options for panel stats are mutually exclusive! Will use the highest ranked choice (panel > panel_intervals > all panels)"
}
if (params.panel) {
	panel = params.genomes[params.assembly].panels[params.panel].intervals
	panels = Channel.fromPath(panel)
} else if (params.panel_intervals) {
	Channel.fromPath(params.panel_intervals)
	.ifEmpty { exit 1; "Could not find the specified gene panel (--panel_intervals)" }
	.set { panels }
} else if (params.all_panels) {
	panel_list = []
	panel_names = params.genomes[params.assembly].panels.keySet()
	panel_names.each {
		interval = params.genomes[params.assembly].panels[it].intervals
		panel_list << file(interval)
	}
	panels = Channel.fromList(panel_list)
} else {
	panels = Channel.empty()
}

// A single exon only covered in male samples - simple sex check
SRY_BED  = params.sry_bed ?: params.genomes[params.assembly].sry_bed
SRY_REGION = file(SRY_BED)

if (!SRY_REGION.exists() ) {
	exit 1, "Could not find the bed file for SRY check!"
}

// Whether to produce BAM output instead of CRAM
params.cram = false
align_suffix = (params.cram == false) ? "bam" : "cram"

// Location of output files
OUTDIR = file(params.outdir)

// Available exome kits

if (TARGETS == false || BAITS == false ) {
   exit 1, "Information on enrichment kit incomplete or missing (please see the documentation for details!)"
}

// Whether to send a notification upon workflow completion
params.email = false

if(params.email == false) {
	exit 1, "You must provide an Email address to which pipeline updates are send!"
}

// Check if the max_length argument is a number
if (params.max_length) {
	summary['maxReadLength'] = params.max_length
 	if (params.max_length instanceof Integer) {

	} else {
		exit 1, "Defined a max_length option, but it is not an integer...!"
	}
}

	
summary['runName'] = run_name
summary['Samples'] = inputFile
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Assembly'] = FASTA
summary['Kit'] = TARGETS
if (params.panel) {
	summary['GenePanel'] = params.panel
} else if (params.panel_intervals) {
	summary['GenePanel'] = params.panel_intervals
} else if (params.all_panels) {
	summary['GenePanel'] = "All panels"
}
summary['CommandLine'] = workflow.commandLine
if (KILL) {
        summary['KillList'] = KILL
}
if (workflow.containerEngine) {
	summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['References'] = [:]
summary['References']['DBSNP'] = DBSNP
summary['IntervallPadding'] = params.interval_padding
summary['SessionID'] = workflow.sessionId

// Header log info
log.info "========================================="
log.info "Exome-seq pipeline v${params.version}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "Exome kit:			${params.kit}"
if (params.panel) {
	log.info "Panel(s):			${params.panel}"
} else if (params.panel_intervals) {
	log.info "Panel(s):			custom"
} else if (params.all_panels) {
	log.info "Panel(s)			all"
} 
log.info "-----------------------------------------"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${run_name}"
if (workflow.containerEngine) {
	log.info "Container engine:		${workflow.containerEngine}"
}
log.info "========================================="

// ******************
// Set up DV channels
// ******************

process list_to_bed {

	input:
        file(targets) from targets_to_bed

        output:
        file(bed) into (bedToExamples, BedToMerge)

        script:
        bed = targets.getBaseName() + ".bed"

        """
                picard IntervalListTools I=$targets O=targets.padded.interval_list PADDING=$params.interval_padding
                picard IntervalListToBed I=targets.padded.interval_list O=$bed
        """
}

faiToExamples = Channel
	.fromPath(FAI_F)
	.ifEmpty{exit 1, "Fai file not found"}

fastaGz = Channel
	.fromPath(FASTAGZ_F)
	.ifEmpty{exit 1, "Fastagz file not found"}
	.into {fastaGzToExamples; fastaGzToVariants }

gzFai = Channel
	.fromPath(GZFAI_F)
	.ifEmpty{exit 1, "gzfai file not found"}
	.into{gzFaiToExamples; gzFaiToVariants }

gzi = Channel
	.fromPath(GZI_F)
	.ifEmpty{exit 1, "gzi file not found"}
	.into {gziToExamples; gziToVariants}

// Read sample file 
Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  readPairsFastp }

process runFastp {

	scratch true

	input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, date, fastqR1, fastqR2 from readPairsFastp

	output:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, date, center, file(left),file(right) into inputBwa
	set file(html),file(json) into fastp_results
	
	script:

	def options = ""
	if (params.max_length != false) {
		options += " -b ${params.max_length} -B ${params.max_length}"
	}

	left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
	json = file(fastqR1).getBaseName() + ".fastp.json"
	html = file(fastqR1).getBaseName() + ".fastp.html"

	"""
		fastp $options --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runBWA {

	scratch true	

	input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa
    
	output:
	set indivID, sampleID, file(outfile) into runBWAOutput
    
	script:
	outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"	
    
	"""
		bwa mem -H $DICT -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${FASTA}\\tCN:${center}" \
			-t ${task.cpus} ${FASTA} $left $right \
			| samtools sort -n -@ 8 -m 3G -O bam -o $outfile - 
	"""	
}

process runFixmate {

	scratch true

	input:
	set indivID, sampleID, file(bam) from runBWAOutput

	output:
	set indivID, sampleID, file(fixed_bam) into FixedBam

	script:
	fixed_bam = bam.getBaseName() + ".fm.bam"

	"""
		samtools fixmate -@2 -m $bam - | samtools sort -@2 -m 3G -O bam -o $fixed_bam -
	"""
}

runBWAOutput_grouped_by_sample = FixedBam.groupTuple(by: [0,1])

process mergeBamFiles_bySample {

	scratch true

	input:
        set indivID, sampleID, file(aligned_bam_list) from runBWAOutput_grouped_by_sample

	output:
	set indivID,sampleID,file(merged_bam),file(merged_bam_index) into mergedBamFile_by_Sample, MergedBamSkipDedup

	script:
	merged_bam = indivID + "_" + sampleID + ".merged.bam"
	merged_bam_index = merged_bam + ".bai"

	if (aligned_bam_list.size() > 1 && aligned_bam_list.size() < 1000 ) {
		"""
			samtools merge -@ 4 $merged_bam ${aligned_bam_list.join(' ')}
			samtools index $merged_bam
		"""
	} else {
		"""
			cp $aligned_bam_list $merged_bam
        	        samtools index $merged_bam
		"""
	}
}

process runMarkDuplicates {

        publishDir "${OUTDIR}/${indivID}/${sampleID}/", mode: 'copy'

      //  scratch true

        input:
        set indivID, sampleID, file(merged_bam),file(merged_bam_index) from mergedBamFile_by_Sample

        output:
        set indivID, sampleID, file(outfile_bam),file(outfile_bai) into BamMD, BamForMultipleMetrics, runHybridCaptureMetrics, runPrintReadsOutput_for_OxoG_Metrics, Bam_for_HC_Metrics, inputPanelCoverage
	set file(outfile_bam), file(outfile_bai) into BamForSexCheck
	file(outfile_md5)
	file(outfile_metrics) into DuplicatesOutput_QC
        val(sample_name) into SampleNames

        script:
        outfile_bam = indivID + "_" + sampleID + ".dedup.bam"
       	outfile_bai = indivID + "_" + sampleID + ".dedup.bam.bai"
	outfile_md5 = indivID + "_" + sampleID + ".dedup.bam.md5"

	sample_name = indivID + "_" + sampleID

        outfile_metrics = indivID + "_" + sampleID + "_duplicate_metrics.txt"

	"""
		samtools markdup -@ ${task.cpus} $merged_bam $outfile_bam
		samtools index $outfile_bam
		samtools stats $outfile_bam > $outfile_metrics
		md5sum $outfile_bam > $outfile_md5
	"""

}

// a simple sex check looking at coverage of the SRY gene
process runSexCheck {

	input:
	file(bams) from BamForSexCheck.collect()

	output:
	file(sex_check_yaml) into SexChecKYaml

	script:
	sex_check_yaml = "sex_check_mqc.yaml"

	"""
		parse_sry_coverage.pl --fasta $FASTA --region $SRY_REGION > $sex_check_yaml
	"""	
}

// ************************
// Run DeepVariant
// ************************

process runDeepvariant {

	label 'deepvariant'

        publishDir "${params.outdir}/${indivID}/${sampleID}/DeepVariant", mode: 'copy'

        input:
        set indivID, sampleID, file(bam),file(bai) from BamMD
        file(bed) from bedToExamples.collect()
        file fai from faiToExamples.collect()
        file fastagz from fastaGzToExamples.collect()
        file gzfai from gzFaiToExamples.collect()
        file gzi from gziToExamples.collect()

        output:
        set indivID,sampleID,file(gvcf)
        file(gvcf) into MergeGVCF
        file(vcf)

        script:
        gvcf = bam.getBaseName() + ".g.vcf.gz"
        vcf = bam.getBaseName() + ".vcf.gz"

        """
                /opt/deepvariant/bin/run_deepvariant \
                --model_type=WES \
                --ref=$fastagz \
                --reads $bam \
                --output_vcf=$vcf \
                --output_gvcf=$gvcf \
                --regions=$bed \
                --num_shards=${task.cpus}
        """
}

if (params.joint_calling) {
	
	process runMergeGvcf {

		label 'glnexus'

                scratch true

                input:
                file(gvcfs) from MergeGVCF.collect()
                file(bed) from BedToMerge.collect()

                output:
                file(merged_vcf) into MergedVCF

                script:
                merged_vcf = "deepvariant.merged." + run_name + ".vcf.gz"

                """
                        /usr/local/bin/glnexus_cli \
                        --config DeepVariantWES \
                        --bed $bed \
                        $gvcfs | bcftools view - | bgzip -c > $merged_vcf

                """
        }

        process annotateIDs {

                publishDir "${OUTDIR}/DeepVariant", mode: 'copy'

                input:
                file (vcf) from MergedVCF

                output:
                set file(vcf_annotated), file(vcf_annotated_index) into VcfAnnotated

                script:
                vcf_annotated = vcf.getBaseName() + ".rsids.vcf.gz"
		vcf_annotated_index = vcf_annotated + ".tbi"

                """
			echo "##reference=${params.assembly}" > header.txt
                        tabix $vcf
                        bcftools annotate -h header.txt -c ID -a $DBSNP -O z -o $vcf_annotated $vcf
			tabix $vcf_annotated
                """
        }

	process VcfGetSample {

		publishDir "${params.outdir}/DeepVariant", mode: 'copy'

		label 'gatk'

                input:
                set file(vcf),file(vcf_index) from VcfAnnotated
                val(sample_name) from SampleNames

                output:
                set file(vcf_sample),file(vcf_sample_index) into VcfSample

                script:
                vcf_sample = sample_name + ".vcf.gz"
		vcf_sample_index = vcf_sample + ".tbi"

                """
			gatk SelectVariants --remove-unused-alternates --exclude-non-variants -V $vcf -sn $sample_name -O $vcf_sample -OVI
                """

        }

        process VcfStats {

                input:
                set file(vcf),file(tbi) from VcfSample

                output:
                file(vcf_stats) into VcfInfo

                script:
                vcf_stats = vcf.getBaseName() + ".stats"

                """
                        bcftools stats $vcf > $vcf_stats
                """

        }

} else {
	VcfInfo = Channel.empty()
}

// *********************
// Compute statistics for fastQ files, libraries and samples
// *********************

process runCollectMultipleMetrics {

	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

	input:
	set indivID, sampleID, bam, bai from BamForMultipleMetrics
        file(baits) from BaitsToMetrics.collect()

	output:
	file("${prefix}*") into CollectMultipleMetricsOutput mode flatten

	script:       
	prefix = indivID + "_" + sampleID + "."

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
		REFERENCE_SEQUENCE=${FASTA} \
		DB_SNP=${DBSNP} \
		INTERVALS=${baits} \
		ASSUME_SORTED=true \
		QUIET=true \
		OUTPUT=${prefix} \
		TMP_DIR=tmp
	"""
}	

process runHybridCaptureMetrics {

	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

	input:
	set indivID, sampleID, file(bam), file(bai) from Bam_for_HC_Metrics
	file(targets) from TargetsToHS.collect()
	file(baits) from BaitsToHS.collect()

	output:
	file(outfile) into HybridCaptureMetricsOutput mode flatten
	file(outfile_per_target)

	script:
	outfile = indivID + "_" + sampleID + ".hybrid_selection_metrics.txt"
	outfile_per_target = indivID + "_" + sampleID + ".hybrid_selection_per_target_metrics.txt"

	"""
        picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
		PER_TARGET_COVERAGE=${outfile_per_target} \
                TARGET_INTERVALS=${targets} \
                BAIT_INTERVALS=${baits} \
                REFERENCE_SEQUENCE=${FASTA} \
		MINIMUM_MAPPING_QUALITY=$params.min_mapq \
                TMP_DIR=tmp
        """
}

process runOxoGMetrics {

    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_OxoG_Metrics
    file(targets) from TargetsToOxo.collect()

    output:
    file(outfile) into runOxoGMetricsOutput mode flatten

    script:
    outfile = indivID + "_" + sampleID + ".OxoG_metrics.txt"

    """

         picard -Xmx${task.memory.toGiga()}G CollectOxoGMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                DB_SNP=${DBSNP} \
                INTERVALS=${targets} \
                REFERENCE_SEQUENCE=${FASTA} \
                TMP_DIR=tmp
        """
}

// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

// this is not finished yet, need to create a proper yaml file
process get_software_versions {

    publishDir "${OUTDIR}/Summary/versions", mode: 'copy'

    output:
    file("v*.txt")
    file(yaml_file) into (software_versions_yaml_fastqc, software_versions_yaml_lib, software_versions_yaml_sample)

    script:
    yaml_file = "software_versions_mqc.yaml"

    """
    echo $workflow.manifest.version &> v_ikmb_exoseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastp -v &> v_fastp.txt
    echo "Deepvariant 1.0.0" &> v_deepvariant.txt
    echo "GLNexus 1.2.6" &> v_glnexus.txt
    samtools --version &> v_samtools.txt
    bcftools --version &> v_bcftools.txt
    multiqc --version &> v_multiqc.txt
    bwa &> v_bwa.txt 2>&1 || true
    parse_versions.pl >  $yaml_file
    """
}

process runMultiqcFastq {

    publishDir "${OUTDIR}/Summary/Fastq", mode: 'copy'

    when:
    !params.skip_multiqc	    

    input:
    file('*') from fastp_results.flatten().toList()
    file('*') from software_versions_yaml_fastqc.collect()
    
    output:
    file("fastp_multiqc*")
    	
    script:

    """
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    cp $params.logo . 
    multiqc -n fastp_multiqc *.json *.html *.yaml
    """
}

process runMultiqcLibrary {

    publishDir "${OUTDIR}/Summary/Library", mode: 'copy'

    when:
    !params.skip_multiqc
	    
    input:
    file('*') from DuplicatesOutput_QC.flatten().toList()
    file('*') from software_versions_yaml_lib.collect()

    output:
    file("library_multiqc*")
    	
    script:

    """
    cp $params.logo .
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n library_multiqc *.txt *.yaml
    """
}

process runMultiqcSample {

    publishDir "${OUTDIR}/Summary/Sample", mode: 'copy'

    when:
    !params.skip_multiqc
	    
    input:
    file('*') from CollectMultipleMetricsOutput.flatten().toList()
    file('*') from HybridCaptureMetricsOutput.flatten().toList()
    file('*') from runOxoGMetricsOutput.flatten().toList()
    file('*') from software_versions_yaml_sample.collect() 
    file('*') from SexChecKYaml.collect()
    file('*') from VcfInfo.collect()
        
    output:
    file("sample_multiqc.html") into multiqc_report
    file("*.yaml")
	
    script:

    def subject = 'Diagnostic exome analysis quality report'
    def recipient = params.email

    """
    cp $params.logo . 
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n sample_multiqc *

    """
}

panel_coverage_data = inputPanelCoverage.combine(panels)

process runPanelCoverage {

	publishDir "${OUTDIR}//Summary/Panel/PanelCoverage", mode: "copy"

        input:
        set indivID,sampleID,file(bam),file(bai),file(panel) from panel_coverage_data

        output:
        set val(panel_name),file(coverage) into outputPanelCoverage
	set indivID,sampleID,file(target_coverage_xls)
	file(target_coverage)

        script:
        panel_name = panel.getSimpleName()
        coverage = indivID + "_" + sampleID + "." +  panel_name  + ".hs_metrics.txt"
	target_coverage = indivID + "_" + sampleID + "." +  panel_name  + ".per_target.hs_metrics.txt"
	target_coverage_xls = indivID + "_" + sampleID + "." + panel_name + ".per_target.hs_metrics_mqc.xlsx"

	// optionally support a kill list of known bad exons
	def options = ""
	if (KILL) {
		options = "--ban ${KILL}"
	}

        // do something here - get coverage and build an XLS sheet
	// First we identify which analysed exons are actually part of the exome kit target definition. 
        """

		picard -Xmx${task.memory.toGiga()}G IntervalListTools \
			INPUT=$panel \
			SECOND_INPUT=$TARGETS \
			ACTION=SUBTRACT \
			OUTPUT=overlaps.interval_list

                picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                        INPUT=${bam} \
                        OUTPUT=${coverage} \
                        TARGET_INTERVALS=${panel} \
                        BAIT_INTERVALS=${panel} \
			CLIP_OVERLAPPING_READS=false \
                        REFERENCE_SEQUENCE=${FASTA} \
                        TMP_DIR=tmp \
	                MINIMUM_MAPPING_QUALITY=$params.min_mapq \
			PER_TARGET_COVERAGE=$target_coverage

		target_coverage2xls.pl $options --infile $target_coverage --min_cov $params.panel_coverage --skip overlaps.interval_list --outfile $target_coverage_xls

        """
}

grouped_panels = outputPanelCoverage.groupTuple()

process runMultiqcPanel {

	publishDir "${OUTDIR}/Summary/Panel", mode: "copy"

	input:
	set val(panel_name),file('*') from grouped_panels

	output:
	file("${panel_name}_multiqc.html")

	script:

	"""
		cp $params.logo . 
		cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
		multiqc -n ${panel_name}_multiqc *
	"""
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['Samples'] = params.samples
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['kit'] = TARGETS
  email_fields['assembly'] = FASTA
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "Diagnostic exome analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[IKMB ExoSeq] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[IKMB ExoSeq] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}

