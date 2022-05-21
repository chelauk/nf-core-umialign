process PICARD_MERGEBAMALIGNMENT {
    tag "$meta.id"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(aligned_unmarked_bam), path(unaligned_marked_bam)
    path(fasta)
    path(dict)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_records = task.memory.toGiga() * 250000
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MergeBamAlignment] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    [ ! -d "./tmpdir" ] && mkdir ./tmpdir || echo "./tmpdir exists"

    picard \\
        -Xmx${avail_mem}g \\
        MergeBamAlignment \\
        TMP_DIR=./tmpdir \\
        MAX_RECORDS_IN_RAM=${max_records} \\
        R=${fasta} \\
        $args \\
        UNMAPPED_BAM=${unaligned_marked_bam} \\
        ALIGNED_BAM=${aligned_unmarked_bam} \\
        ADD_MATE_CIGAR=true CLIP_ADAPTERS=false \\
        CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true \\
        MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\
        ATTRIBUTES_TO_RETAIN=XS \\
        OUTPUT=${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$( echo \$(picard MergeBamAlignment --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}" 
    """
    touch ${prefix}_pre_collapse.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 2.26.10 
    END_VERSIONS
    """

}
