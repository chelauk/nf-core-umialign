process FGBIO_GROUPREADSBYUMI {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fgbio=2.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.0--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(taggedbam)

    output:
    tuple val(meta), path("*_umi-grouped.bam")  , emit: bam
    tuple val(meta), path("*_umi_histogram.txt"), emit: histogram
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -d "./tmpdir" ] && mkdir ./tmpdir || echo "./tmpdir exists"

    
    fgbio \\
        -Xmx${task.memory.toGiga()}g \\
        -XX:+AggressiveOpts -XX:+AggressiveHeap \\
        --tmp-dir=./tmpdir \\
        GroupReadsByUmi \\
        $args \\
        -i $taggedbam \\
        -o ${prefix}_umi-grouped.bam \\
        -f ${prefix}_umi_histogram.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args > args.txt
    touch ${prefix}_umi-grouped.bam
    touch ${prefix}_umi_histogram.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 1.3 
    END_VERSIONS
    """
}
