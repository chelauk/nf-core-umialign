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
                meta.id = meta.patient + "-" + meta.sample
                [[meta.patient, meta.sample, meta.id, meta.gender, meta.status],bam ]}
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
                                    [meta,bam]
                                }
        // STEP 1.5: MERGING AND INDEXING BAM FROM MULTIPLE LANES
        SAMBAMBA_MERGE(bam_multiple)
        bam       = bam_single.mix(SAMBAMBA_MERGE.out.bam)

    emit:
        bam
}
