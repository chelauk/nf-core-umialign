process FGBIO_ERROR_RATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.3.0--0' :
        'quay.io/biocontainers/fgbio:1.3.0--0' }"

    input:
    tuple val(meta), path(bam)
    path interval_list
    path fasta
    path dict
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path("*txt")           , emit: error_rate
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgbio -Xmx${task.memory.toGiga()}g \\
        ErrorRateByReadPosition \\
        --input $bam \\
        --variants $dbsnp \\
        --intervals $interval_list \\
        --ref $fasta \\
        --output ${prefix}_${args} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}_${args}_error_rate.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: 1.3 
    END_VERSIONS
    """
}
