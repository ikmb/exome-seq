#!/usr/bin/env nextflow

/**
===============================
Exome Pipeline
===============================

This Pipeline performs variant calling
using Google DeepVariant

### Homepage / git
git@github.com:ikmb/exome-seq.git
### Implementation
Re-Implemented in Q3 2020

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
--cnv 			       Enable calling of copy number variants (for select assembly and kit combinations)
--phase			       Perform phasing of the final call set. 
--joint_calling		       Perform joint calling of all samples (default: true)
--amplicon		       This is a small amplicon-based analysis, skip duplicate marking and sex check
--skip_multiqc		       Don't attached MultiQC report to the email. 
--panel 		       Gene panel to check coverage of (valid options: cardio_dilatative, cardio_hypertrophic, cardio_non_compaction, eoIBD_25kb, imm_eoIBD_full, breast_cancer)
--all_panels 		       Run all gene panels defined for this assembly (none if no panel is defined!)
--panel_intervals	       Run a custom gene panel in interval list format (must have a matching sequence dictionary!)
--run_name 		       A descriptive name for this pipeline run
--cram			       Whether to output the alignments in CRAM format (default: bam)
--interval_padding	       Include this number of nt upstream and downstream around the exome targets (default: 10)
--vep			       Perform variant annotation with VEP (requires substantial local configuration work!)
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
--cnv_ref		       A CNVkit reference  (ccn.gz) that matches the assembly and kit used. Use with care. 
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
BWA2_INDEX = file(params.genomes[ params.assembly ].bwa2_index)

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

if (TARGETS==BAITS) {
	exit 1, "Target and bait files must not have the same name to avoid file collisions!"
}
targets_to_bed = Channel.fromPath(TARGETS)


Channel.from(file(TARGETS))
	.into { TargetsToHS; TargetsToMetrics; TargetsToOxo; TargetsToPanel }

Channel.from(file(BAITS))
	.into { BaitsToHS; BaitsToMetrics }

if (params.kill) {
	KILL = params.kill
} else if (params.kit && params.genomes[params.assembly].kits[params.kit].kill) {
	KILL = params.genomes[params.assembly].kits[params.kit].kill
} else {
	KILL = false
}

// CNVkit reference
if (params.cnv) {

	if (params.cnv_ref) {

		cnv_ref_file = file(params.cnv_ref).getName()

		Channel.fromPath(params.cnv_ref)
			.ifEmpty { exit 1; "Could not find the specified CNV reference" }
			.set { cnv_ref_gz }

	} else if ( params.genomes[params.assembly].kits[params.kit].containsKey("cnvkit") ) {
		cnv_ref_file = params.genomes[params.assembly].kits[params.kit].cnvkit
		Channel.fromPath(cnv_ref_file)
			.ifEmpty { exit 1; "Could not find a CNVkit reference for this kit and assembly" }
			.set { cnv_ref_gz }

	} else {
		exit 1, "Requested to run CNVkit but no CNV reference is defined for this assembly and exome kit."
	}

} else {
	cnv_ref_gz = Channel.empty()
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

if (params.amplicon) {
	log.info "This is an amplicon run - no duplicate marking or sex-checking will be performed!"
}

if (params.vep) {

	if (params.assembly != "GRCh38") {
		exit 1, "VEP is not currently set up to work with defunct assembly versions...please use GRCh38"
	}

	if (!params.vep_cache_dir || !params.vep_plugin_dir) {
		exit 1, "Missing VEP cache and/or plugin directory..."
	}

	if (params.dbnsfp_db) {
		dbNSFP_DB = file(params.dbnsfp_db)
		if (!dbNSFP_DB.exists()) {
			exit 1, "Could not find the specified dbNSFP database..."
		}
	} else {
		exit 1, "No dbNSFP database defined for this execution profile..."
	}

	if (params.dbscsnv_db) {
		dbscSNV_DB = file(params.dbscsnv_db)
		if ( !dbscSNV_DB.exists() ) {
			exit 1, "Could not find the specified dbscSNV database..."
		}
	} else {
		exit 1, "No dbscSNV database defined for this execution profile..."
	}

	if (params.cadd_snps && params.cadd_indels) {
		CADD_SNPS = file(params.cadd_snps)
		CADD_INDELS = file(params.cadd_indels)
		if (!CADD_SNPS.exists() || !CADD_INDELS.exists() )
			exit 1, "Missing CADD SNPs and/or Indel references..."
		}
	else {
		exit 1, "CADD SNP and/or Indel reference files not defined for this execution profile..."
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
summary['Phasing'] = params.phase
summary['AmpliconRun'] = params.amplicon
summary['CommandLine'] = workflow.commandLine
if (KILL) {
        summary['KillList'] = KILL
}
if (params.cnv_ref || params.cnv) {
	summary["CNVkit CNN"] = cnv_ref_file
}
if (workflow.containerEngine) {
	summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['References'] = [:]
summary['References']['DBSNP'] = DBSNP
if (params.vep) {
	summary['References']['dbNSFP'] = params.dbnsfp_db
	summary['References']['dbSCSNV'] = params.dbscsnv_db
	summary['References']['CADD_SNPs'] = params.cadd_snps
	summary['References']['CADD_Indels'] = params.cadd_indels
}
summary['IntervallPadding'] = params.interval_padding
summary['SessionID'] = workflow.sessionId

// Header log info
log.info "========================================="
log.info "Exome-seq pipeline v${params.version}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "Exome kit:			${params.kit}"
if (params.amplicon) {
	log.info "Amplicon run:			true"
}
if (params.panel) {
	log.info "Panel(s):			${params.panel}"
} else if (params.panel_intervals) {
	log.info "Panel(s):			custom"
} else if (params.all_panels) {
	log.info "Panel(s)			all"
}
if (params.vep) {
	log.info "Run VEP				${params.vep}"
} 
log.info "CNVCalling			${params.cnv}"
log.info "Phasing				${params.phase}"
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

// Deepvariant needs a BED file, convert
process list_to_bed {

	executor 'local' 

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

// Read trimming
process trim {

	scratch params.scratch

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
		fastp -c $options --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

// alignment
process align {

	//scratch true	

	input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa
    
	output:
	set indivID, sampleID, file(outfile) into FixedBam
    
	script:
	outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.fm.bam"	
    
	"""
		bwa-mem2 mem -K 1000000 -H $DICT -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${FASTA}\\tCN:${center}" \
			-t ${task.cpus} ${BWA2_INDEX} $left $right \
			| samtools fixmate -@ ${task.cpus} -m - - \
			| samtools sort -@ ${task.cpus} -m 3G -O bam -o $outfile - 
	"""	
}

// Samples may run on one or more lanes; combine based on keys
FixedBam
	.groupTuple(by: [0,1])
	.into { bams_for_merging;  bams_singleton }

// If sample has multiple bam files, merge
process merge_multi_lane {

	scratch params.scratch

	input:
        set indivID, sampleID, file(aligned_bam_list) from bams_for_merging.filter { i,s,b -> b.size() > 1 && b.size() < 1000 }

	output:
	set indivID,sampleID,file(merged_bam),file(merged_bam_index) into merged_bams

	script:
	merged_bam = indivID + "_" + sampleID + ".merged.bam"
	merged_bam_index = merged_bam + ".bai"
        sample_name = indivID + "_" + sampleID

	"""
			samtools merge -@ 4 $merged_bam ${aligned_bam_list.join(' ')}
	"""
}

// combine merged bam files with singleton bams 
all_bams = merged_bams.concat(bams_singleton.filter { i,s,b -> b.size() < 2 || b.size() > 1000 } )

// index the lot
process bam_index {

	scratch params.scratch

	publishDir "${OUTDIR}/${indivID}/${sampleID}/", mode: 'copy', enabled: params.amplicon

	input:
        set indivID, sampleID, file(bam) from all_bams

	output:
	set indivID, sampleID, file(bam) into bam_indexed

	script:
	bam_index = bam.getName() + ".bai"

	"""
		samtools index $bam
	"""

}

// if this is amplicon data, don't do deduping
if (params.amplicon) {

	bam_indexed
	.into { BamMD; BamForMultipleMetrics; runHybridCaptureMetrics; runPrintReadsOutput_for_OxoG_Metrics; Bam_for_HC_Metrics; inputPanelCoverage ; Bam_for_Cnv; BamForSexCheck ; BamPhasing}

	DuplicatesOutput_QC = Channel.from(false)

	SexChecKYaml = Channel.from(false)
	
} else {

	process dedup {

	        publishDir "${OUTDIR}/${indivID}/${sampleID}/", mode: 'copy'

	      //  scratch true

        	input:
	        set indivID, sampleID, file(merged_bam),file(merged_bam_index) from bam_indexed

        	output:
	        set indivID, sampleID, file(outfile_bam),file(outfile_bai) into BamMD, BamForMultipleMetrics, runHybridCaptureMetrics, runPrintReadsOutput_for_OxoG_Metrics, Bam_for_HC_Metrics, inputPanelCoverage, Bam_for_Cnv
		set file(outfile_bam), file(outfile_bai) into (BamForSexCheck, BamPhasing )
		file(outfile_md5)
		file(outfile_metrics) into DuplicatesOutput_QC

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
	process sex_check {

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

}

// ************************
// Run DeepVariant
// ************************

process deepvariant {

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
        set indivID,sampleID,file(vcf) into Vcf_to_Cnv
	val(sample_name) into SampleNames

        script:
        gvcf = bam.getBaseName() + ".g.vcf.gz"
        vcf = bam.getBaseName() + ".vcf.gz"
	sample_name = indivID + "_" + sampleID

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
	
	process merge_gvcf {

		label 'glnexus'

                scratch params.scratch

                input:
                file(gvcfs) from MergeGVCF.collect()
                file(bed) from BedToMerge.collect()

                output:
                file(merged_vcf) into MergedVCF

                script:
                merged_vcf = "deepvariant.merged." + run_name + ".vcf.gz"

                """
                        /usr/local/bin/glnexus_cli \
                        --config ${params.glnexus_config} \
                        --bed $bed \
                        $gvcfs | bcftools view - | bgzip -c > $merged_vcf

                """
        }

	if (params.phase) {
		process whatshap {
		
			label 'whatshap'

			publishDir "${OUTDIR}/DeepVariant/Phased", mode: 'copy'

			input:
			file (vcf) from MergedVCF
			file('*') from BamPhasing.collect()

			output:
			file(phased_vcf) into PhasedVcf

			script:
			phased_vcf = vcf.getBaseName() + ".phased.vcf.gz"

			"""
				whatshap phase -o $phased_vcf --tag=PS --reference $FASTA $vcf *.bam
			"""
		

		}

	} else {
		PhasedVcf = MergedVCF
	}

        process annotateIDs {

                publishDir "${OUTDIR}/DeepVariant", mode: 'copy'

                input:
                file (vcf) from PhasedVcf

                output:
                set file(vcf_annotated), file(vcf_annotated_index) into VcfAnnotated, VcfToVep

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

	process vep {

		label 'vep'

		publishDir "${OUTDIR}/DeepVariant/VEP", mode: 'copy'

		when:
		params.vep

		input:
		set file(vcf),file(vcf_index) from VcfToVep

		output:
		file(vcf_vep)
		file(vcf_alissa)

		script:
		vcf_vep = vcf.getBaseName() + ".vep.vcf"
		vcf_alissa = vcf.getBaseName() + ".vep2alissa.vcf"

		"""
			vep --offline \
				--cache \
				--dir ${params.vep_cache_dir} \
				--species homo_sapiens \
				--assembly $params.assembly \
				-i $vcf \
				--format vcf \
				-o $vcf_vep --dir_plugins ${params.vep_plugin_dir} \
				--plugin dbNSFP,$dbNSFP_DB,${params.dbnsfp_fields} \
				--plugin dbscSNV,$dbscSNV_DB \
				--plugin CADD,${params.cadd_snps},${params.cadd_indels} \
				--plugin ExACpLI \
				--fasta $FASTA \
				--fork 4 \
				--vcf \
				--per_gene \
				--sift p \
				--polyphen p \
				--check_existing \
				--canonical
	
			sed -i.bak 's/CADD_PHRED/CADD_phred/g' $vcf_vep

			vep2alissa.pl --infile $vcf_vep > $vcf_alissa
		"""
	}

	process vcf_get_sample {

		//publishDir "${params.outdir}/DeepVariant", mode: 'copy'

		label 'gatk'

                input:
                set file(vcf),file(vcf_index) from VcfAnnotated
                val(sample_name) from SampleNames

                output:
                set file(vcf_sample),file(vcf_sample_index) into VcfSample, VcfReheader

                script:
                vcf_sample = sample_name + ".vcf.gz"
		vcf_sample_index = vcf_sample + ".tbi"

                """
			gatk SelectVariants --remove-unused-alternates --exclude-non-variants -V $vcf -sn $sample_name -O variants.vcf.gz -OVI
		        gatk LeftAlignAndTrimVariants -R $FASTA -V variants.vcf.gz -O $vcf_sample
			rm variants.vcf.gz

                """

        }

	process vcf_add_header {

                publishDir "${params.outdir}/DeepVariant", mode: 'copy'

		input:
		set file(vcf),file(tbi) from VcfReheader

		output:
		set file(vcf_r),file(tbi_r) 

		script:

		vcf_r = vcf.getBaseName() + ".final.vcf.gz"
		tbi_r = vcf_r + ".tbi"

		"""
			echo "##reference=${params.assembly}" > header.txt
			bcftools annotate -h header.txt -O z -o $vcf_r $vcf
			tabix $vcf_r
		"""

	}

        process vcf_stats {

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


// *******************************************
// Optional: CNVkit with pre-defined reference
// *******************************************
if (params.cnv) {

        process stage_cnv_reference {

                executor 'local'

                input:
                file(ref_gz) from cnv_ref_gz

                output:
                file(ref_cnn) into cnv_ref

                script:
                ref_cnn = ref_gz.getBaseName()

                """
                        gunzip -c $ref_gz > $ref_cnn
                """

        }

        process cnvkit_batch{

                label 'cnvkit'

                publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Processing", mode: 'copy'

                input:
                set indivID, sampleID, file(bam),file(bai) from Bam_for_Cnv
                file(cnn) from cnv_ref.collect()

                output:
                set indivID, sampleID,file(cnr),file(cns) into Cnv_to_seg

                script:
                cnr = bam.getBaseName() + ".cnr"
                cns = bam.getBaseName() + ".cns"
                """
                        cnvkit.py batch -r $cnn -d out *.bam
                        mv out/*.cns .
                        mv out/*.cnr .
                """

        }

        process cnvkit_segmetrics {
                label 'cnvkit'

                publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Processing", mode: 'copy'

                input:
                set indivID, sampleID, file(cnr),file(cns) from Cnv_to_seg

                output:
                set indivID, sampleID,file(cnr),file(seg_cns) into Cnv_to_call

                script:
                seg_cns = cns.getBaseName() + ".segmetrics.cns"

                """
                        cnvkit.py segmetrics -s $cns $cnr --ci
                """

        }

	Cnv_call_vcf = Cnv_to_call.join(Vcf_to_Cnv, by: [0,1] )

        process cnvkit_call {

                label 'cnvkit'

                publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Processing", mode: 'copy'

                input:
                set indivID, sampleID, file(cnr),file(cns),file(vcf) from Cnv_call_vcf

                output:
                set indivID, sampleID,file(cnr),file(call_cns) into Cnv_to_gene, Cnv_to_break, Cnv_to_export, Cnv_to_plot

                script:
                call_cns = cns.getBaseName() + ".call.cns"

                """
                        cnvkit.py call $cns --filter ci
                """

        }

        process cnvkit_genemetrics {

                label 'cnvkit'

                publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Metrics", mode: 'copy'

                input:
                set indivID, sampleID, file(cnr),file(cns) from Cnv_to_gene

                output:
                file(metrics)

                script:

                metrics = cnr.getBaseName() + ".genemetrics.txt"

                """
                        cnvkit.py genemetrics -s $cns $cnr -t 0.2 > $metrics
                """

        }

        process cnvkit_breaks {

                label 'cnvkit'

                publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Metrics", mode: 'copy'

                input:
                set indivID, sampleID, file(cnr),file(cns) from Cnv_to_break

                output:
                file(breaks)

                script:

                breaks = cnr.getBaseName() + ".breaks.txt"

                """
                        cnvkit.py breaks $cns $cnr > $breaks
                """

        }

	process cnvkit_export {
	
		label 'cnvkit'

		publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit", mode: 'copy'

		input:
		set indivID, sampleID,file(cnr),file(call_cns) from Cnv_to_export

		output:
		set file(bed),file(vcf) into CnvOut

		script:
		bed = call_cns.getBaseName() + ".bed"
		vcf = call_cns.getBaseName() + ".vcf"

		"""
			cnvkit.py export bed $call_cns -o $bed
			cnvkit.py export vcf $call_cns -i $sampleID -o $vcf
		"""

	}

	process cnvkit_plots {
	
		label 'cnvkit'

		publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Plots", mode: 'copy'

		input:
                set indivID, sampleID,file(cnr),file(call_cns) from Cnv_to_plot

		output:
		file(scatter)
		file(diagram)

		script:
		scatter = call_cns.getBaseName() + ".scatter.pdf"
		diagram = call_cns.getBaseName() + ".diagram.pdf"

		"""
			cnvkit.py scatter --y-min -4 --y-max 4 -o $scatter -s $call_cns $cnr
			cnvkit.py diagram -o $diagram -s $call_cns $cnr 
		"""
	}
}

// *********************
// Compute statistics for fastQ files, libraries and samples
// *********************

process multi_metrics {

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

process hybrid_capture_metrics {

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

process oxo_metrics {

    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_OxoG_Metrics
    file(targets) from TargetsToOxo.collect()

    output:
    file(outfile) into runOxoGMetricsOutput mode flatten

    script:
    outfile = indivID + "_" + sampleID + ".OxoG_metrics.txt"

    """

         picard -Xmx${task.memory.toGiga()-1}G CollectOxoGMetrics \
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
    echo "CNVkit 0.9.9" &> v_cnvkit.txt
    echo "Deepvariant 1.2.0" &> v_deepvariant.txt
    echo "GLNexus 1.3.1" &> v_glnexus.txt
    samtools --version &> v_samtools.txt
    bcftools --version &> v_bcftools.txt
    multiqc --version &> v_multiqc.txt
    bwa-mem2 &> v_bwa.txt 2>&1 || true
    parse_versions.pl >  $yaml_file
    """
}

process multiqc_fastq {

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

process multiqc_library {

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
    multiqc -n library_multiqc *
    """
}

process multiqc_sample {

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

process panel_coverage {

	publishDir "${OUTDIR}//Summary/Panel/PanelCoverage", mode: "copy"

        input:
        set indivID,sampleID,file(bam),file(bai),file(panel) from panel_coverage_data
	file(targets) from TargetsToPanel.collect()

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
			SECOND_INPUT=$targets \
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

process multiqc_panel {

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

