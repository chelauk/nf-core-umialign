//
// stage one:
// fastqs converted to bam and split into two channels
// the first channel is aligned to the reference
// then merged with the unaligned bam to include the umi info
// fastqc, library complexity, adapter contamination and error rate
// metrics are also obtained.
//

include { FASTQC                            } from '../../../modules/nf-core/modules/fastqc/main'
include { MULTIQC                           } from '../../../modules/nf-core/modules/multiqc/main'
include { FGBIO_FASTQTOBAM                  } from '../../../modules/nf-core/modules/fgbio/fastqtobam/main'
include { PICARD_ESTIMATELIBRARYCOMPLEXITY  } from '../../../modules/nf-core/modules/picard/estimatelibrarycomplexity/main'
include { PICARD_MARKADAPTERS               } from '../../../modules/nf-core/modules/picard/markadapters/main'
include { PICARD_BAMTOFASTQ as B2FQ1        } from '../../../modules/nf-core/modules/picard/bamtofastq/main'
include { BWA_MEM as BM1                    } from '../../../modules/nf-core/modules/bwa/mem/main'
include { PICARD_SORTSAM as ALIGNED_SORT1   } from '../../../modules/nf-core/modules/picard/sortsam/main' 
include { PICARD_SORTSAM as UNALIGNED_SORT1 } from '../../../modules/nf-core/modules/picard/sortsam/main' 
include { PICARD_MERGEBAMALIGNMENT as PMB1  } from '../../../modules/nf-core/modules/picard/mergebamalignment/main'
include { PICARD_COLLECTHSMETRICS as HS1    } from '../../../modules/nf-core/modules/picard/collecthsmetrics/main'
include { FGBIO_ERROR_RATE as ER1           } from '../../../modules/nf-core/modules/fgbio/errorrate/main'


// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules

include { MERGE_RUNS                        } from '../merge_runs/main'

workflow UMI_STAGE_ONE {
    take:
    intervals
    reads      // channel: [ val(meta)]
    fasta
    dict
    dbsnp
    dbsnp_tbi
    bwa

    main:
    ch_versions = Channel.empty()

    // MODULE: Run FastQC

    FASTQC ( reads )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // MODULE: fastq to bam

    FGBIO_FASTQTOBAM ( reads, params.read_structure )
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())

    // MODULE: library complexity

    PICARD_ESTIMATELIBRARYCOMPLEXITY ( FGBIO_FASTQTOBAM.out.umibam )
    ch_versions = ch_versions.mix(PICARD_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())

    // MODULE: Picard mark illumina adapters

    PICARD_MARKADAPTERS (FGBIO_FASTQTOBAM.out.umibam)
    ch_versions = ch_versions.mix(PICARD_MARKADAPTERS.out.versions.first())

    // MODULE: Picard bam to fastq

    B2FQ1 (PICARD_MARKADAPTERS.out.bam)
    ch_versions = ch_versions.mix(B2FQ1.out.versions.first())

    // MODULE: bwa align to reference

    sort = true // unnecessary but helps to conform to nf-core template
    BM1 ( B2FQ1.out.fastq,bwa,sort )
    ch_versions = ch_versions.mix(BM1.out.versions.first())

    // MODULE: sort before merge
    aligned_suffix= "aligned_namesorted"
    ALIGNED_SORT1(BM1.out.bam, aligned_suffix)

    // MODULE: sort before merge
    unaligned_suffix= "unaligned_namesorted"
    UNALIGNED_SORT1(PICARD_MARKADAPTERS.out.bam,unaligned_suffix)
    
    // MODULE: merge adapter marked unaligned bam to aligned bam
    PMB1 ( ALIGNED_SORT1.out.bam.join(UNALIGNED_SORT1.out.bam),fasta,dict )
    ch_versions = ch_versions.mix(PMB1.out.versions.first())

    // SUBWORKFLOW: merge bams if necessary

    MERGE_RUNS(PMB1.out.bam)

    // MODULE : Picard collect hs metrics before collapse

    HS1 ( MERGE_RUNS.out.bam, intervals )
    ch_versions = ch_versions.mix(HS1.out.versions.first())

    // MODULE : FGBIO collect error rates prior to collapse

    ER1 ( MERGE_RUNS.out.bam, intervals, fasta, dict, dbsnp, dbsnp_tbi )
    ch_versions = ch_versions.mix(ER1.out.versions.first())

    emit:
    complexity_metrics   = PICARD_ESTIMATELIBRARYCOMPLEXITY.out.complexity_metrics // channel: [ val(meta), txt ]
    adapter_metrics      = PICARD_MARKADAPTERS.out.metrics                         // channel: [ val(meta), txt ]
    merged_bam           = MERGE_RUNS.out.bam                                      // channel: [ val(meta), bam ]
    fastqc_log           = FASTQC.out.zip                                          // channel: [ val(meta), zip ]
    pre_collapse_metrics = HS1.out.hs_metrics                                      // channel: [ val(meta), txt ]
    pre_collapse_error   = ER1.out.error_rate                                      // channel: [ val(meta), txt ]

    versions           = ch_versions                                               // channel: [ versions.yml ]

}
