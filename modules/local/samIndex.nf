#!/usr/bin/env nextflow

process samIndex {

    input:
        path bamPath

    output:
        path "${bamPath}.bai"
    
    script:
    """
    samtools index $bamPath
    """
}