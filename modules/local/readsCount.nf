#!/usr/bin/env nextflow

process readsCount {

    input: 
        path readsPath

    output:
        path "${readsPath.baseName}_read_count.csv"
    
    script:
    """
    echo "${readsPath.baseName}, \$(samtools import $readsPath | samtools view -c)" >> ${readsPath.baseName}_read_count.csv
    """    
}