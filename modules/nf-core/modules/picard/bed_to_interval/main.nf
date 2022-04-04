process PICARD_BED_TO_INTERVAL_LIST {
    tag "interval_list"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    path target_bed
    path dict

    output:
    path "interval.list"                          , emit: intervals
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard ] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        BedToIntervalList \\
        --INPUT $target_bed \\
        --SEQUENCE_DICTIONARY $dict \\
        --OUTPUT interval.list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard BedToIntervalList --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    """
    touch interval.list
    touch versions.yml
    """
}
