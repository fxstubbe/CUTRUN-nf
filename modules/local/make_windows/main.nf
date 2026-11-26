process MAKE_WINDOWS {
    tag "${wmeta.bin}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedtools:2.31.1--a120a7e98287539a' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

    input:
    tuple val(wmeta), path(fai)     // wmeta: [bin: '10kb', size: 10000]
    val exclude_chroms              // comma-separated chromosome names to leave out ('' = none)

    output:
    tuple val(wmeta), path("*.windows.bed"), emit: bed
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${params.genome}.${wmeta.bin}"
    """
    # Windows follow the order of the FASTA index, which is also the BAM header order,
    # so bedtools coverage can run with -sorted. The last window of each chromosome is
    # truncated at the chromosome end.
    awk -v OFS='\\t' -v ex='${exclude_chroms}' \\
        'BEGIN { n = split(ex, a, ","); for (i = 1; i <= n; i++) drop[a[i]] = 1 }
         !(\$1 in drop) { print \$1, \$2 }' ${fai} > genome.sizes

    bedtools makewindows -g genome.sizes -w ${wmeta.size} > ${prefix}.windows.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
