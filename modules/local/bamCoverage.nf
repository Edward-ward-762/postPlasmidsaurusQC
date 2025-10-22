#!/usr/bin/env nextflow

params.bamPath='bamPath_local_placeholder'

process bamCoverage {
    
    publishDir 'output', mode: 'copy'

    input: 
        path bamPath

    output:
        path "${bamPath.baseName}_coverage.bedgraph"
    
    script:
    """
    samtools index $bamPath
    bamCoverage -b $bamPath -of bedgraph -bs 10000 -o "${bamPath.baseName}_coverage.bedgraph"
    """
}