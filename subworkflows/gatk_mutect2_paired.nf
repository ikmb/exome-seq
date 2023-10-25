include { GATK_MUTECT2_PAIR }                   from "./../modules/gatk/mutect2_pair"
include { GATK_SPLITINTERVALS }                 from './../modules/gatk/splitintervals'
include { GATK_FILTER_MUTECT_CALLS }            from "./../modules/gatk/filter_mutect_calls"
include { BCFTOOLS_VIEW }                       from "./../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP }             from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_ANNOTATE }                   from "./../modules/bcftools/annotate"
include { GATK_LEARN_READ_ORIENTATION_MODEL }   from "./../modules/gatk/learn_read_orientation_model"
include { GATK_GET_PILEUP_SUMMARIES }           from "./../modules/gatk/get_pileup_summaries"
include { GATK_CALCULATE_CONTAMINATION_PAIRED } from "./../modules/gatk/calculate_contamination_paired"
include { GATK_MERGEMUTECTSTATS }               from "./../modules/gatk/mergemutectstats"
include { GATK_MERGEVCFS }                      from "./../modules/gatk/mergevcfs"

ch_versions     = Channel.from([])
ch_vcfs         = Channel.from([])

workflow GATK_MUTECT2_PAIRED {

    take:
        bams
        targets
        fasta
        dbsnp
        mutect_normals
        mutect_normals_tbi

    main:

        // Make split targets
        GATK_SPLITINTERVALS(
            targets,
            fasta
        )
        GATK_SPLITINTERVALS.out.intervals.flatMap { i ->
            i.collect { file(it) }
        }.set { targets_split }

        // extract all the normal bams from the per-patient data drop
        ch_normal_bam = bams.map { m,nb,nbi,tb,tbi ->
            [[
            patient_id: m.patient_id,
            sample_id: m.normal_id,
            original_sample_id: m.sample_id,
            status: 0
            ], nb, nbi ]
        }

        // get all the tumor bams from the per-patient data drop
        ch_tumor_bam = bams.map { m,nb,nbi,tb,tbi ->
            [ m,tb,tbi ]
        }.transpose()

        ch_tumor_bam_clean = ch_tumor_bam.map { m,tbam,tbi ->
            [[
            patient_id: m.patient_id,
            sample_id: tbam.getName().split("-dedup")[0],
            status: 1
            ],tbam,tbi ]
        }

        // combine all bams into a serial emission of individual bam files
        ch_all_bams = ch_normal_bam.mix(ch_tumor_bam_clean)

        // Produce raw mutect2 multi-sample vcf against normal(s)
        // This call is parallelized by calling chunk (split intervals)
        GATK_MUTECT2_PAIR(
            bams.collect(),
            targets_split,
            fasta,
            mutect_normals,
            mutect_normals_tbi
        )

        ch_versions = ch_versions.mix(GATK_MUTECT2_PAIR.out.versions)

        // Merge all chunked stats into one stat file per patient
        GATK_MERGEMUTECTSTATS(
            GATK_MUTECT2_PAIR.out.stats.groupTuple()
        )

        ch_versions = ch_versions.mix(GATK_MERGEMUTECTSTATS.out.versions)

        // Merge Mutect2 calls across chunks
        // The groupTuple produces [ meta, [ vcfs ], [ tbis ]]
        GATK_MERGEVCFS(
            GATK_MUTECT2_PAIR.out.vcf.map{ m,v,t,s -> [m,v,t] }.groupTuple()
        )

        ch_versions = ch_versions.mix(GATK_MERGEVCFS.out.versions)
        
        // Learn the read orientation model across all chunks
        // This is grouped across samples
        GATK_LEARN_READ_ORIENTATION_MODEL(
            GATK_MUTECT2_PAIR.out.f1r2.groupTuple()
        )

        ch_versions = ch_versions.mix(GATK_LEARN_READ_ORIENTATION_MODEL.out.versions)
        
        // Get pileup summaries for each BAM file and each target (not chunk!)
        GATK_GET_PILEUP_SUMMARIES(
            ch_all_bams,
            targets,
            fasta
        )
        
        ch_versions = ch_versions.mix(GATK_GET_PILEUP_SUMMARIES.out.versions)

        // Split pileup summaries by tumor status
        GATK_GET_PILEUP_SUMMARIES.out.table.branch { m,t ->
            normal: m.status == 0
            tumor: m.status = 1
        }.set { ch_pileup_status }
        
        ch_versions = ch_versions.mix(GATK_GET_PILEUP_SUMMARIES.out.versions)

        ch_pileup_normal = ch_pileup_status.normal.map { m,t ->
            [ m.patient_id,m,t ]
        }

        ch_pileup_tumor = ch_pileup_status.tumor.map { m,t ->
            [ m.patient_id,m,t ]
        }

        // combine pipleup summaries for normal samples with each respective tumor sample
        ch_pileup_joined = ch_pileup_normal.cross(ch_pileup_tumor)
        
        ch_pileup_joined.map { normal,tumor ->
            [[
				patient_id: normal[1].patient_id,
				normal_id: normal[1].sample_id,
				tumor_id: tumor[1].sample_id,
				sample_id: "${tumor[1].sample_id}_vs_${normal[1].sample_id}",
                original_sample_id: normal[1].original_sample_id
			],normal[2],tumor[2]]
        }.set { ch_pileup_pairs }
        
        // calculate contamination for each tumor-normal pairing
        GATK_CALCULATE_CONTAMINATION_PAIRED(
            ch_pileup_pairs
        )
        
        ch_versions = ch_versions.mix(GATK_CALCULATE_CONTAMINATION_PAIRED.out.versions)

        ch_contamtination_grouped = GATK_CALCULATE_CONTAMINATION_PAIRED.out.table.map { m,t ->
            [[
                patient_id: m.patient_id,
                sample_id: m.original_sample_id,
                normal_id: m.normal_id,
                tumor_id: m.tumor_id
            ],t]
        }

        ch_versions = ch_versions.mix(GATK_CALCULATE_CONTAMINATION_PAIRED.out.versions)
        
        // merge all results together into a per-patient data package again (incl. grouping 1-n contamination tables)
	
        GATK_MERGEVCFS.out.vcf.map { m,v,t ->
            [ m.patient_id,m,v,t ]
        }
        .join(
            GATK_MERGEMUTECTSTATS.out.stats.map { m,s ->
                [ m.patient_id,s ]
            }    
        ).join(
            GATK_LEARN_READ_ORIENTATION_MODEL.out.model.map { m,f -> 
                [ m.patient_id, f ]
            }      
        ).join(
            ch_contamtination_grouped.map { m,t ->
                [ m.patient_id, t]
            }.groupTuple()
        ).map { k,m,v,t,s,f,tbl ->
            [ m,v,t,s,f,tbl ]
        }.set { ch_mutect }
        
        // Finally filter the original multi-sample VCF using the contamination and strand information produced above
        GATK_FILTER_MUTECT_CALLS(
            ch_mutect,
            fasta
        )

        ch_versions = ch_versions.mix(GATK_FILTER_MUTECT_CALLS.out.versions)

        ch_vcfs = ch_vcfs.mix(GATK_FILTER_MUTECT_CALLS.out.vcf)

        BCFTOOLS_VIEW(ch_vcfs)

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        BCFTOOLS_ANNOTATE_DBSNP(
            BCFTOOLS_VIEW.out.vcf.map { meta,v,t ->
                def s_meta = [ id: meta.id, sample_id: meta.sample_id, patient_id: meta.patient_id, variantcaller: "MUTECT2" ]
                tuple(s_meta,v,t)
            },
            dbsnp
        )
        
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_DBSNP.out.versions)

        BCFTOOLS_ANNOTATE(
            BCFTOOLS_ANNOTATE_DBSNP.out.vcf
        )

    emit:
    versions     = ch_versions
    vcf         = BCFTOOLS_ANNOTATE.out.vcf
        
}

