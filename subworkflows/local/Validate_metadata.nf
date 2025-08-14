
process VALIDATE_METADATA_PROCESS {
    input:
    path metadata_csv
    path validate_script
    output:
    path 'validation_report.txt', emit: validation_report
    script:
    """
    if [ ! -f "$metadata_csv" ]; then
        echo "ERROR: Metadata CSV file does not exist: $metadata_csv" >&2
        exit 1
    fi
    python3 $validate_script $metadata_csv > validation_report.txt
    """
}

workflow VALIDATE_METADATA {
    take:
    metadata_csv

    main:
    script_ch = Channel.value('../modules/local/validate_metadata.py')
    validation = VALIDATE_METADATA_PROCESS(metadata_csv, script_ch)

    emit:
    validation
}
