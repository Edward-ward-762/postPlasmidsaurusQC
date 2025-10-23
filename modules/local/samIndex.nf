#!/usr/bin/env nextflow

process samIndex {

    input:
        tuple val(meta), path(bamPath)

    output:
        tuple val(meta), path("${bamPath}.bai"), emit: bai
    
    script:
    """
    samtools index $bamPath
    """
}