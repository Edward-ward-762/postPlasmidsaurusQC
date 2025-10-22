#!/usr/bin/env nextflow

params.readsPath='reads_local_placeholder'

process readsCount {
    
    publishDir 'output', mode: 'copy'

    input: 
        path readsPath

    output:
        path "${readsPath.baseName}_read_count.csv"
    
    script:
    """
    echo "${readsPath.baseName}, \$(samtools import $readsPath | samtools view -c)" >> ${readsPath.baseName}_read_count.csv
    """    
}