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

if (params.input) csv_file = file(params.input)
ch_input_sample = extract_csv(csv_file)
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
include { MERGE_RUNS  } from '../subworkflows/local/merge_runs/main'
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
    ch_bam_st2  = Channel.empty()

    //
    // MODULE: Target bed to interval list
    //

    PICARD_BED_TO_INTERVAL_LIST (
        target_bed,dict
    )

    if (params.stage == 'one') {
        ch_input_sample.branch{
            fastq: it[0].data_type == "fastq"
            bam:   it[0].data_type == "bam"
        }.set{ch_input_sample_type}

        //
        // MODULE: Run FastQC
        //

        FASTQC (
            ch_input_sample_type.fastq
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        //
        // MODULE: fastq to bam
        //

        FGBIO_FASTQTOBAM (
            ch_input_sample_type.fastq,
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
        ch_versions = ch_versions.mix(PMB1.out.versions.first())

        //
        // SUBWORKFLOW: merge bams if necessary
        //

        MERGE_RUNS(PMB1.out.bam)

        //
        // MODULE : Picard collect hs metrics before collapse
        //

        HS1 (
            MERGE_RUNS.out.bam,
            PICARD_BED_TO_INTERVAL_LIST.out.intervals
        )

        //
        // MODULE : FGBIO collect error rates prior to collapse
        //

        ER1 (
            MERGE_RUNS.out.bam,
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
            MERGE_RUNS.out.bam
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


        if (params.stage == "two" ) {
            ch_bam_st2 = ch_input_sample
        }
            ch_bam_st2 = ch_input_sample_type.bam.mix(FGBIO_FILTERCONSENSUSREADS.out.bam)

            //
            // MODULE: Picard bam to fastq
            //

            B2FQ2 (
                ch_bam_st2
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
                BM2.out.bam.join(SAMBAMBA_SORT.out.bam),
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
            // MODULE : FGBIO collect error rates after collapse
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
            ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

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
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]

        //TODO since it is mandatory: error/warning if not present?
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no gender specified, gender is not considered
        // gender is only mandatory for somatic CNV
        if (row.gender) meta.gender = row.gender.toString()
        else meta.gender = "NA"

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        // mapping with fastq
        if (row.lane && row.fastq_3) {
            meta.patient    = row.patient.toString()
            meta.sample     = row.sample.toString()
            meta.lane       = row.lane.toString()
            meta.id         = "${row.patient}_${row.sample}_${row.lane}"
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def fastq_3     = file(row.fastq_3, checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.patient}_${row.sample}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.patient}_${row.sample}\\tPL:ILLUMINA\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'
            return [meta, [fastq_1, fastq_2, fastq_3]]
        // start from BAM
        } else if (row.lane && row.bam) {
            meta.id         = "${row.sample}_${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.patient}_${row.sample}\\tPL:ILLUMINA\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = "bam"
            return [meta, bam]
        // recalibration
        } else if (row.table && row.cram) {
            meta.id   = meta.sample
            def cram  = file(row.cram,  checkIfExists: true)
            def crai  = file(row.crai,  checkIfExists: true)
            def table = file(row.table, checkIfExists: true)
            meta.data_type  = "cram"
            return [meta, cram, crai, table]
        // recalibration when skipping MarkDuplicates
        } else if (row.table && row.bam) {
            meta.id   = meta.sample
            def bam   = file(row.bam,   checkIfExists: true)
            def bai   = file(row.bai,   checkIfExists: true)
            def table = file(row.table, checkIfExists: true)
            meta.data_type  = "bam"
            return [meta, bam, bai, table]
        // prepare_recalibration or variant_calling
        } else if (row.cram) {
            meta.id = meta.sample
            def cram = file(row.cram, checkIfExists: true)
            def crai = file(row.crai, checkIfExists: true)
            meta.data_type  = "cram"
            return [meta, cram, crai]
        // prepare_recalibration when skipping MarkDuplicates
        } else if (row.bam) {
            meta.id = meta.sample
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)
            meta.data_type  = "bam"
            return [meta, bam, bai]
        // annotation
        } else if (row.vcf) {
            meta.id = meta.sample
            def vcf = file(row.vcf, checkIfExists: true)
            meta.data_type  = "vcf"
            return [meta, vcf]
        } else {
            log.warn "Missing or unknown field in csv file header"
        }
    }
}

// Function to count number of intervals
def count_intervals(intervals_file) {
    count = 0

    intervals_file.eachLine{ it ->
        count += it.startsWith("@") ? 0 : 1
    }

    return count
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
