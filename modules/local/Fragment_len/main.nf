process FRAGMENT_LEN {
    tag "$meta.id"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt") , emit: tsv
    path  "versions.yml"           , emit: versions

    when:
    params.fragment_len_reporting  // only run if true

    script:
    def args     = task.ext.args ?: ''
    def args2    = task.ext.args2 ?: ''
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    samtools view $args -@ $task.cpus $bam | $args2 > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
