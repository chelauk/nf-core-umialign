//
// stage three:
// fgbio filter, return bam to fastq realign with bwa
// get post collapse metrics
//

include { FGBIO_FILTERCONSENSUSREADS        } from '../../../modules/nf-core/modules/fgbio/filterconsensus/main'
include { PICARD_BAMTOFASTQ as B2FQ2        } from '../../../modules/nf-core/modules/picard/bamtofastq/main'
include { BWA_MEM as BM2                    } from '../../../modules/nf-core/modules/bwa/mem/main'
include { PICARD_MERGEBAMALIGNMENT as PMB2  } from '../../../modules/nf-core/modules/picard/mergebamalignment/main'
include { PICARD_COLLECTHSMETRICS as HS2    } from '../../../modules/nf-core/modules/picard/collecthsmetrics/main'
include { FGBIO_ERROR_RATE as ER2           } from '../../../modules/nf-core/modules/fgbio/errorrate/main'
include { PICARD_UMIMARKDUPLICATES          } from '../../../modules/nf-core/modules/picard/umimarkduplicates/main'
include { PICARD_COLLECTINSERTSIZEMETRICS   } from '../../../modules/nf-core/modules/picard/collectinsertsizemetrics/main'
include { QUALIMAP_BAMQC                    } from '../../../modules/nf-core/modules/qualimap/bamqc/main'

workflow UMI_STAGE_THREE {

    take:
    intervals
    bam
    fasta
    bwa
    dict
    target_bed
    dbsnp
    dbsnp_tbi

    main:
    ch_versions = Channel.empty()

    // MODULE: filter consensus reads

    FGBIO_FILTERCONSENSUSREADS ( bam,fasta )
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())

    // MODULE: Picard bam to fastq

    B2FQ2 ( FGBIO_FILTERCONSENSUSREADS.out.bam )
    ch_versions = ch_versions.mix( B2FQ2.out.versions.first())

    // MODULE: second bwa alignment
    sort = false
    BM2 ( B2FQ2.out.fastq,bwa,sort )

    // MODULE: picard mergebamalignment
    PMB2(BM2.out.bam.join(FGBIO_FILTERCONSENSUSREADS.out.bam),fasta,dict)
    ch_versions = ch_versions.mix(PMB2.out.versions.first())

    // MODULE : Picard collect hs metrics after collapse

    HS2 ( PMB2.out.bam, intervals )
    ch_versions = ch_versions.mix(HS2.out.versions.first())

    // MODULE : FGBIO collect error rates after collapse

    ER2 ( PMB2.out.bam,intervals,fasta,dict,dbsnp,dbsnp_tbi )
    ch_versions = ch_versions.mix(ER2.out.versions.first())

    // MODULE: PICARD UMI-AWARE markduplicates

    PICARD_UMIMARKDUPLICATES ( PMB2.out.bam )
    ch_versions = ch_versions.mix(PICARD_UMIMARKDUPLICATES.out.versions.first())
    
   // MODULE: PICARD collect insert size metrics

    PICARD_COLLECTINSERTSIZEMETRICS  ( PMB2.out.bam )
    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())

    // MODULE: Qualimap BamQC

    QUALIMAP_BAMQC ( PICARD_UMIMARKDUPLICATES.out.bam,target_bed )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    emit:
    md_umi_metrics        = PICARD_UMIMARKDUPLICATES.out.umi_metrics          // channel: [ val(meta), txt ]
    md_metrics            = PICARD_UMIMARKDUPLICATES.out.md_metrics           // channel: [ val(meta), txt ]
    post_collapse_metrics = HS2.out.hs_metrics                                // channel: [ val(meta), txt ]
    post_collapse_error   = ER2.out.error_rate                                // channel: [ val(meta), txt ]
    bamqc                 = QUALIMAP_BAMQC.out.results                        // channel: [ val(meta), folder ]
    insert_sizes          = PICARD_COLLECTINSERTSIZEMETRICS.out.size_metrics  // channel: [ val(meta), txt ]
    versions              = ch_versions                                       // channel: [ versions.yml ]
}
