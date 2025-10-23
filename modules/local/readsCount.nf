#!/usr/bin/env nextflow

process readsCount {
    label 'process_low'

    input: 
        path readsPath

    output:
        path "${readsPath.baseName}_read_count.csv", emit: count
        path  "versions.yml"                       , emit: versions
    
    script:
    """
    echo "${readsPath.baseName}, \$(samtools import $readsPath | samtools view -c)" >> ${readsPath.baseName}_read_count.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1)
    END_VERSIONS
    """    
}