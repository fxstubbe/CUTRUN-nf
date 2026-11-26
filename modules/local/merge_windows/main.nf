process MERGE_WINDOW_COUNTS {
    tag "${bin}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(bin), path(coverage_files, stageAs: 'coverage/*')
    path spikein_logs, stageAs: 'spikein_log/*'
    path metadata

    output:
    tuple val(bin), path("window_counts_${bin}.csv.gz"), emit: counts
    tuple val(bin), path("library_metrics_${bin}.csv") , emit: metrics
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    merge_window_counts.py \\
        --coverage-dir coverage \\
        --spikein-log-dir spikein_log \\
        --metadata ${metadata} \\
        --normalisation-c ${params.normalisation_c} \\
        --out-counts window_counts_${bin}.csv.gz \\
        --out-metrics library_metrics_${bin}.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
