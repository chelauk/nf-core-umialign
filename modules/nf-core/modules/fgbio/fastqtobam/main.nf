process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fgbio=1.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.4.0--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:1.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    val read_structure

    output:
    tuple val(meta), path("*_unaln.bam"), emit: umibam
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir tmp

    fgbio \\
        -Xmx${task.memory.toGiga()}g \\
        --tmp-dir=${PWD}/tmp \\
        FastqToBam \\
        $args \\
        -i $reads \\
        --sort true \\
        -o "${prefix}_unaln.bam" \\
        --read-structures $read_structure \\
        --read-group-id $meta.lane \\
        --umi-tag RX \\
        --sample $meta.id \\
        --library $meta.id 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unaln.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 1.3
    END_VERSIONS
    """
}
