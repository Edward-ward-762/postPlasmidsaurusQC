#!/usr/bin/env nextflow

process bamCoverage {

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