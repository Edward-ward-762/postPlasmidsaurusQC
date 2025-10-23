#!/usr/bin/env nextflow

process bamCoverage {

    input: 
        tuple val(meta), path(bamPath), path(baiPath)

    output:
        tuple val(meta), path("${bamPath.baseName}_coverage.bedgraph"), emit: bedgraph
    
    script:
    """
    bamCoverage -b $bamPath -of bedgraph -bs $params.bedgraph_bin_size -o "${bamPath.baseName}_coverage.bedgraph"
    """
}