#!/usr/bin/env nextflow

process collectReadCounts {
    label 'process_single'

    input:
        path input_files

    output:
        path "read_count.csv"

    script:
    """
    cat ${input_files} > "read_count.csv"
    """
}