#!/usr/bin/env nextflow

params.bamPath='bamPath_local_placeholder'

process nanoPlot {
    
    publishDir "output/${readsPath.baseName}", mode: 'copy'

    container 'tbc'
    
    input: 
        path readsPath

    output:
        path "${readsPath.baseName}_NanoPlot"
    
    script:
    """
    NanoPlot --fastq $readsPath -o ${readsPath.baseName}_NanoPlot
    """
}