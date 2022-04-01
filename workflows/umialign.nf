/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowUmialign.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [  params.input,
                            params.multiqc_config,
                            params.fasta,
                            params.bwa
                        ]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// initialise refs and indices indices
fasta       = params.fasta                      ? Channel.fromPath(params.fasta).collect() : Channel.empty()
bwa         = params.fasta ? params.bwa         ? Channel.fromPath(params.bwa).collect()   : null : []
dict        = params.fasta ? params.dict        ? Channel.fromPath(params.dict).collect()  : null : []
dbsnp       = params.fasta ? params.dbsnp       ? Channel.fromPath(params.dbsnp).collect() : null : []
dbsnp_tbi   = params.fasta ? params.dbsnp_tbi   ? Channel.fromPath(params.dbsnp_tbi).collect() : null : []

// initialise params outside the genome scope
target_bed  = params.target_bed    ? file(params.target_bed) : file("${params.outdir}/no_file")
file("${params.outdir}/no_file").text = "no_file\n"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                            } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main'
include { CAT_FASTQ                         } from '../modules/nf-core/modules/cat/fastq/main'
include { PICARD_BED_TO_INTERVAL_LIST       } from '../modules/nf-core/modules/picard/bed_to_interval/main'
include { FGBIO_FASTQTOBAM                  } from '../modules/nf-core/modules/fgbio/fastqtobam/main'
include { PICARD_ESTIMATELIBRARYCOMPLEXITY  } from '../modules/nf-core/modules/picard/estimatelibrarycomplexity/main'
include { PICARD_MARKADAPTERS               } from '../modules/nf-core/modules/picard/markadapters/main'
include { PICARD_BAMTOFASTQ as B2FQ1        } from '../modules/nf-core/modules/picard/bamtofastq/main'
include { BWA_MEM as BM1                    } from '../modules/nf-core/modules/bwa/mem/main'
include { PICARD_MERGEBAMALIGNMENT as PMB1  } from '../modules/nf-core/modules/picard/mergebamalignment/main'
include { PICARD_COLLECTHSMETRICS as HS1    } from '../modules/nf-core/modules/picard/collecthsmetrics/main'
include { FGBIO_ERROR_RATE as ER1           } from '../modules/nf-core/modules/fgbio/errorrate/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { FGBIO_GROUPREADSBYUMI             } from '../modules/nf-core/modules/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../modules/nf-core/modules/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        } from '../modules/nf-core/modules/fgbio/filterconsensus/main'
include { PICARD_BAMTOFASTQ as B2FQ2        } from '../modules/nf-core/modules/picard/bamtofastq/main'
include { BWA_MEM as BM2                    } from '../modules/nf-core/modules/bwa/mem/main'
include { SAMBAMBA_SORT                     } from '../modules/nf-core/modules/sambamba/sort/main'
include { PICARD_MERGEBAMALIGNMENT as PMB2  } from '../modules/nf-core/modules/picard/mergebamalignment/main'
include { PICARD_COLLECTHSMETRICS as HS2    } from '../modules/nf-core/modules/picard/collecthsmetrics/main'
include { FGBIO_ERROR_RATE as ER2           } from '../modules/nf-core/modules/fgbio/errorrate/main'
include { PICARD_UMIMARKDUPLICATES          } from '../modules/nf-core/modules/picard/umimarkduplicates/main'
include { QUALIMAP_BAMQC                    } from '../modules/nf-core/modules/qualimap/bamqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow UMIALIGN {

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map { meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // MODULE: Target bed to interval list
    //

    PICARD_BED_TO_INTERVAL_LIST (
        target_bed,dict
    )

    //
    // MODULE: Run FastQC
    //

    FASTQC (
        ch_cat_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: fastq to bam
    //

    FGBIO_FASTQTOBAM (
        ch_cat_fastq,
        params.read_structure
    )
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())

    //
    // MODULE: library complexity
    //

    PICARD_ESTIMATELIBRARYCOMPLEXITY (
        FGBIO_FASTQTOBAM.out.umibam
    )
    ch_versions = ch_versions.mix(PICARD_ESTIMATELIBRARYCOMPLEXITY.out.versions.first())

    //
    // MODULE: Picard mark illumina adapters
    //

    PICARD_MARKADAPTERS (
        FGBIO_FASTQTOBAM.out.umibam
    )
    ch_versions = ch_versions.mix(PICARD_MARKADAPTERS.out.versions.first())

    //
    // MODULE: Picard bam to fastq
    //

    B2FQ1 (
        PICARD_MARKADAPTERS.out.bam
    )
    ch_versions = ch_versions.mix(B2FQ1.out.versions.first())

    //
    // MODULE: bwa align to reference
    //

    sort = true
    BM1 (
        B2FQ1.out.fastq,
        bwa,
        sort
    )
    ch_versions = ch_versions.mix(BM1.out.versions.first())

    //
    // MODULE: merge adapter marked unaligned bam to aligned bam
    //

    PMB1 (
        BM1.out.bam.join(PICARD_MARKADAPTERS.out.bam),
        fasta,
        dict
    )

    //
    // MODULE : Picard collect hs metrics before collapse
    //

    HS1 (
        PMB1.out.bam,
        PICARD_BED_TO_INTERVAL_LIST.out.intervals
    )

    //
    // MODULE : FGBIO collect errorates prior to collapse
    //

    ER1 (
        PMB1.out.bam,
        PICARD_BED_TO_INTERVAL_LIST.out.intervals,
        fasta,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    //
    // MODULE: FGBIO group reads by UMI
    //
    FGBIO_GROUPREADSBYUMI (
        PMB1.out.bam
    )

    //
    // MODULE: FGBIO call consensus
    //

    FGBIO_CALLMOLECULARCONSENSUSREADS (
        FGBIO_GROUPREADSBYUMI.out.bam
    )

    //
    // MODULE: FGBIO filter consensus reads
    //

    FGBIO_FILTERCONSENSUSREADS (
        FGBIO_CALLMOLECULARCONSENSUSREADS.out.bam,fasta
    )

    //
    // MODULE: Picard bam to fastq
    //

    B2FQ2 (
        FGBIO_FILTERCONSENSUSREADS.out.bam
    )

    //
    // MODULE: second alignment
    //

    BM2 (
        B2FQ2.out.fastq,
        bwa,
        sort
    )

    SAMBAMBA_SORT (
        FGBIO_FILTERCONSENSUSREADS.out.bam
    )
    ch_versions = ch_versions.mix(SAMBAMBA_SORT.out.versions.first())

    PMB2 (
        BM1.out.bam.join(SAMBAMBA_SORT.out.bam),
        fasta,
        dict
    )

    //
    // MODULE : Picard collect hs metrics after collapse
    //

    HS2 (
        PMB2.out.bam,
        PICARD_BED_TO_INTERVAL_LIST.out.intervals
    )

    //
    // MODULE : FGBIO collect errorates after collapse
    //

    ER2 (
        PMB2.out.bam,
        PICARD_BED_TO_INTERVAL_LIST.out.intervals,
        fasta,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    PICARD_UMIMARKDUPLICATES (
        PMB2.out.bam
    )

    QUALIMAP_BAMQC (
        PMB2.out.bam,
        target_bed
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowUmialign.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ER1.out.error_rate.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ER2.out.error_rate.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(HS1.out.hs_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(HS2.out.hs_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_UMIMARKDUPLICATES.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_ESTIMATELIBRARYCOMPLEXITY.out.complexity_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKADAPTERS.out.metrics.collect{it[1]}.ifEmpty([]))
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
