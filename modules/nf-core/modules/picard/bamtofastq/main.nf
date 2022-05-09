process PICARD_BAMTOFASTQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.fastq.gz")   , emit: fastq
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def max_records = task.memory.toGiga() * 100000
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    [ ! -d "./tmpdir" ] && mkdir ./tmpdir || echo "./tmpdir exists"

    picard \\
        -Xmx${avail_mem}g \\
        SamToFastq \\
        MAX_RECORDS_IN_RAM=${max_records} \\
        TMP_DIR=./tmpdir \\
        $args \\
        I=$bam \\
        FASTQ=${prefix}.fastq.gz \\
        CLIPPING_ATTRIBUTE=XT \\
        CLIPPING_ACTION=2 \\
        INTERLEAVE=true \\
        NON_PF=true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.26.10 
    END_VERSIONS
    """
}
