#!/usr/bin/env nextflow

params.bamPath='bamPath_local_placeholder'

process samIndex {
    
    publishDir 'output', mode: 'copy'

    input:
        path bamPath

    output:
        path "${bamPath}.bai"
    
    script:
    """
    samtools index $bamPath
    """
}