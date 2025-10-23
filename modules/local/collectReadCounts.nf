#!/usr/bin/env nextflow

process collectReadCounts {
    
    input:
        path input_files

    output:
        path "read_count.csv"

    script:
    """
    cat ${input_files} > "read_count.csv"
    """
}