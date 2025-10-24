#!/usr/bin/env nextflow

process readsCount {
    label 'process_low'

    input: 
        path readsPath

    output:
        path "${readsPath.baseName}_read_count.csv", emit: count
    
    script:
    """
    echo "${readsPath.baseName}, \$(samtools import $readsPath | samtools view -c)" >> ${readsPath.baseName}_read_count.csv
    """    
}