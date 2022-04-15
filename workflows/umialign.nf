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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF_CORE modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PICARD_BED_TO_INTERVAL_LIST       } from '../modules/nf-core/modules/picard/bed_to_interval/main'
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UMI_STAGE_ONE                     } from '../subworkflows/local/umi_stages/umi_st_one'
include { UMI_STAGE_TWO                     } from '../subworkflows/local/umi_stages/umi_st_two'
include { UMI_STAGE_THREE                   } from '../subworkflows/local/umi_stages/umi_st_three'


include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow UMIALIGN {
// initialise channels
    ch_versions  = Channel.empty()
    bam_input1   = Channel.empty()
    bam_input2   = Channel.empty()
	umi_st_one_adapter_metrics         = Channel.empty()
	umi_st_one_complexity_metrics      = Channel.empty()
	umi_st_one_fastqc_log              = Channel.empty()
	umi_st_one_pre_collapse_error      = Channel.empty()
	umi_st_one_pre_collapse_metrics    = Channel.empty()
	umi_st_two_family_sizes            = Channel.empty()
	umi_st_three_md_umi_metrics        = Channel.empty()
	umi_st_three_md_metrics            = Channel.empty()
	umi_st_three_post_collapse_error   = Channel.empty()
	umi_st_three_post_collapse_metrics = Channel.empty()
	umi_st_three_bamqc                 = Channel.empty()


    // MODULE create interval list with picard

    PICARD_BED_TO_INTERVAL_LIST (target_bed,dict)

    // SUBWORKFLOW Run from fastqs

    if ( params.stage == 'one' ) {
        reads = extract_csv(csv_file)
    	UMI_STAGE_ONE( PICARD_BED_TO_INTERVAL_LIST.out.intervals, reads,fasta,dict,dbsnp,dbsnp_tbi,bwa )
    	bam_input1 = bam_input1.mix(UMI_STAGE_ONE.out.merged_bam)
   		umi_st_one_adapter_metrics       = UMI_STAGE_ONE.out.adapter_metrics.mix(UMI_STAGE_ONE.out.adapter_metrics)
    	umi_st_one_complexity_metrics    = UMI_STAGE_ONE.out.complexity_metrics.mix(UMI_STAGE_ONE.out.complexity_metrics)
    	umi_st_one_fastqc_log            = UMI_STAGE_ONE.out.fastqc_log.mix(UMI_STAGE_ONE.out.fastqc_log)
    	umi_st_one_pre_collapse_error    = UMI_STAGE_ONE.out.pre_collapse_error.mix(UMI_STAGE_ONE.out.pre_collapse_error)
    	umi_st_one_pre_collapse_metrics  = UMI_STAGE_ONE.out.pre_collapse_metrics.mix(UMI_STAGE_ONE.out.pre_collapse_metrics)
		ch_versions = ch_versions.mix(UMI_STAGE_ONE.out.versions.first())
    }

    // MODULE: Input merged bams before group and call consensus

    if ( params.stage == 'two' ) {
        bam_input1 = extract_csv(csv_file)
        UMI_STAGE_TWO(bam_input1)
        bam_input2  = bam_input2.mix(UMI_STAGE_TWO.out.consensus_bam)
        umi_st_two_family_sizes            = UMI_STAGE_TWO.out.family_sizes.mix(UMI_STAGE_TWO.out.family_sizes)
        ch_versions = ch_versions.mix(UMI_STAGE_TWO.out.versions.first())
    } else {
        UMI_STAGE_TWO(bam_input1)
        bam_input2  = bam_input2.mix(UMI_STAGE_TWO.out.consensus_bam)
        umi_st_two_family_sizes            = UMI_STAGE_TWO.out.family_sizes.mix(UMI_STAGE_TWO.out.family_sizes)
        ch_versions = ch_versions.mix(UMI_STAGE_TWO.out.versions.first())
    }

    // MODULE: Input output of call consensus

    if ( params.stage == 'three' ) {
        bam_input2  = extract_csv(csv_file)
        UMI_STAGE_THREE(PICARD_BED_TO_INTERVAL_LIST.out.intervals,bam_input2,fasta,bwa,dict,target_bed,dbsnp,dbsnp_tbi)
        umi_st_three_md_umi_metrics        =  UMI_STAGE_THREE.out.md_umi_metrics.mix(UMI_STAGE_THREE.out.md_umi_metrics)
        umi_st_three_md_metrics            =  UMI_STAGE_THREE.out.md_metrics.mix(UMI_STAGE_THREE.out.md_metrics)
        umi_st_three_post_collapse_error   =  UMI_STAGE_THREE.out.post_collapse_error.mix(UMI_STAGE_THREE.out.post_collapse_error)
        umi_st_three_post_collapse_metrics =  UMI_STAGE_THREE.out.post_collapse_metrics.mix(UMI_STAGE_THREE.out.post_collapse_metrics)
        umi_st_three_bamqc                 =  UMI_STAGE_THREE.out.bamqc.mix(UMI_STAGE_THREE.out.bamqc)
        ch_versions = ch_versions.mix(UMI_STAGE_THREE.out.versions.first())
    } else {
        UMI_STAGE_THREE(PICARD_BED_TO_INTERVAL_LIST.out.intervals,bam_input2,fasta,bwa,dict,target_bed,dbsnp,dbsnp_tbi)
        umi_st_three_md_umi_metrics        =  UMI_STAGE_THREE.out.md_umi_metrics.mix(UMI_STAGE_THREE.out.md_umi_metrics)
        umi_st_three_md_metrics            =  UMI_STAGE_THREE.out.md_metrics.mix(UMI_STAGE_THREE.out.md_metrics)
        umi_st_three_post_collapse_error   =  UMI_STAGE_THREE.out.post_collapse_error.mix(UMI_STAGE_THREE.out.post_collapse_error)
        umi_st_three_post_collapse_metrics =  UMI_STAGE_THREE.out.post_collapse_metrics.mix(UMI_STAGE_THREE.out.post_collapse_metrics)
        umi_st_three_bamqc                 =  UMI_STAGE_THREE.out.bamqc.mix(UMI_STAGE_THREE.out.bamqc)
        ch_versions = ch_versions.mix(UMI_STAGE_THREE.out.versions.first())
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))

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
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_one_adapter_metrics.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_one_complexity_metrics.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_one_fastqc_log.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_one_pre_collapse_error.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_one_pre_collapse_metrics.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_two_family_sizes.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_three_md_umi_metrics.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_three_md_metrics.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_three_post_collapse_error.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_three_post_collapse_metrics.unique().collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(umi_st_three_bamqc.unique().collect{it[1]}.ifEmpty([]))
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
        if (row.lane) {
            meta.patient    = row.patient.toString()
            meta.sample     = row.sample.toString()
            meta.lane       = row.lane.toString()
            meta.id         = "${row.patient}_${row.sample}_${row.lane}"
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.patient}_${row.sample}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.patient}_${row.sample}\\tPL:ILLUMINA\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'

            fqs = [] // list to fill with all arguements matching regex below
            map_fqs = row.findAll{k,v -> k.matches(~/^fastq(_[0-9]+)?/)} // find fqs in row
            map_fqs = map_fqs.sort{it -> (it.key.replaceAll(/fastq(_)?/, "") ?: 0).toInteger() } // sort matches
            for (fq in map_fqs) {
            fqs.add(file(fq.value, checkIfExists: true))
            }

            return [meta, fqs]

        // start from BAM
        } else if (row.patient && row.bam) {
            meta.patient    = row.patient.toString()
            meta.sample     = row.sample.toString()
            meta.id         = "${row.patient}_${row.sample}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def read_group  = "\"@RG\\tID:${row.patient}_${row.sample}\\tLB:${row.patient}_${row.sample}\\tSM:${row.patient}_${row.sample}\\tPL:ILLUMINA\""
            meta.read_group = read_group.toString()
            meta.data_type  = "bam"
            return [meta, bam]
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
