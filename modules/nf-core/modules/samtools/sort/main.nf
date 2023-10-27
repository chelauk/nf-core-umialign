process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda 'modules/nf-core/samtools/sort/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(aligned_unmarked_bam), path(unaligned_marked_bam)

    output:
    tuple val(meta), path("*aligned_namesorted.bam"), path("*marked_namesorted.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        $args \\
        -n \\
        -@ $task.cpus \\
        -o ${prefix}.aligned_namesorted.bam \\
        -T $prefix \\
        $aligned_unmarked_bam

    samtools sort \\
        $args \\
        -n \\
        -@ $task.cpus \\
        -o ${prefix}.marked_namesorted.bam \\
        -T $prefix \\
        $unaligned_unmarked_bam
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aligned_namesorted.bam
    touch ${prefix}.marked_namesorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
