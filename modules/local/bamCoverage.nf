#!/usr/bin/env nextflow

process bamCoverage {

    input: 
        path bamPath

    output:
        path "${bamPath.baseName}_coverage.bedgraph"
    
    script:
    """
    samtools index $bamPath
    bamCoverage -b $bamPath -of bedgraph -bs $params.bedgraph_bin_size -o "${bamPath.baseName}_coverage.bedgraph"
    """
}