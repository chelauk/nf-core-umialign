//
// stage two:
// fgbio groups umi reads and calls consensus reads
// callconsensus is the slowest stage by far and needs optimisation
// family size metrics are also recorded
//

include { FGBIO_GROUPREADSBYUMI             } from '../../../modules/nf-core/modules/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../../../modules/nf-core/modules/fgbio/callmolecularconsensusreads/main'

workflow UMI_STAGE_TWO {

    take:
    bam

    main:
    ch_versions = Channel.empty()

    // MODULE: FGBIO group reads by UMI

    FGBIO_GROUPREADSBYUMI ( bam )
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())

    // MODULE: FGBIO call consensus

    FGBIO_CALLMOLECULARCONSENSUSREADS ( FGBIO_GROUPREADSBYUMI.out.bam )
    ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions.first())

    emit:
    family_sizes         = FGBIO_GROUPREADSBYUMI.out.histogram            // channel: [ val(meta), txt ]
    consensus_bam        = FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam      // channel: [ val(meta), bam ]

    versions           = ch_versions                                      // channel: [ versions.yml ]
}
