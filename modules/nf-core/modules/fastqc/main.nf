process FASTQC {
    tag "fastqc"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.patient}_${meta.sample}_${meta.lane}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    [ ! -f  ${prefix}_3.fastq.gz ] && ln -s ${reads[2]} ${prefix}_3.fastq.gz
    fastqc $args --threads $task.cpus ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix}_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    // Add soft-links for consistent naming
    def prefix = task.ext.prefix ?: "${meta.patient}_${meta.sample}_${meta.lane}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    [ ! -f  ${prefix}_3.fastq.gz ] && ln -s ${reads[2]} ${prefix}_3.fastq.gz
    touch ${prefix}_1.html
    touch ${prefix}_2.html
    touch ${prefix}_3.html
    touch ${prefix}_1.zip
    touch ${prefix}_2.zip
    touch ${prefix}_3.zip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.11.9
    END_VERSIONS
    """
}
