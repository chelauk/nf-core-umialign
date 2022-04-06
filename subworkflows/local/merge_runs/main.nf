/*
================================================================================
                                MERGING
================================================================================
*/

include { SAMBAMBA_MERGE }              from '../../../modules/nf-core/modules/sambamba/merge/main'

workflow MERGE_RUNS {
    take:
        take: bam_bwa

    main:
    //    bam_bwa = MERGE_BAM_ALIGNMENT.out.bam

        bam_bwa
            .map{meta, bam ->
                CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
                meta.id = meta.patient + "_" + meta.sample
                meta.read_group  = "\"@RG\\tID:1\\t${CN}PU:1\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.patient}_${meta.sample}\\tPL:ILLUMINA\""
                [[meta.patient, meta.sample, meta.id, meta.gender, meta.status, meta.read_group],bam ]}
            .groupTuple()
            .branch{
                single:   it[1].size() == 1
                multiple: it[1].size() > 1
            }.set{ bam_bwa_to_sort }

        bam_multiple = bam_bwa_to_sort.multiple
                                .map{
                                    info,bam ->
                                    def meta = [:]
                                    meta.patient = info[0]
                                    meta.sample  = info[1]
                                    meta.id      = info[2]
                                    meta.gender  = info[3]
                                    meta.status  = info[4]
                                    meta.read_group = info[5]
                                    [meta,bam]
                                }

        bam_single = bam_bwa_to_sort.single
                                .map{
                                    info,bam ->
                                    def meta = [:]
                                    meta.patient = info[0]
                                    meta.sample  = info[1]
                                    meta.id      = info[2]
                                    meta.gender  = info[3]
                                    meta.status  = info[4]
                                    meta.read_group = info[5]
                                    [meta,bam]
                                }
        // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES
        SAMBAMBA_MERGE(bam_multiple)
        bam       = bam_single.mix(SAMBAMBA_MERGE.out.bam)

    emit:
        bam
}
