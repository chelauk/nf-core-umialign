process FGBIO_CALLMOLECULARCONSENSUSREADS {
    tag "$meta.id"
    label 'process_long_medium'

    conda (params.enable_conda ? "bioconda::fgbio=2.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.0--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    [ ! -d "./tmpdir" ] && mkdir ./tmpdir || echo "./tmpdir exists"

    fgbio -Xmx${task.memory.toGiga()}g \\
        -XX:+AggressiveOpts -XX:+AggressiveHeap \\
        --tmp-dir=./tmpdir \\
        CallMolecularConsensusReads \\
        -i $bam \\
        --min-input-base-quality 30 \\
        --read-group-id ${meta.id} \\
        --tag MI \\
        $args \\
        -o ${prefix}.bam
    
    rm -rf./tmpdir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 1.3 
    END_VERSIONS
    """

}
