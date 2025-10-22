#!/usr/bin/env nextflow

params.bamPath='bamPath_local_placeholder'

process fastQC {
    
    publishDir "output/${readsPath.baseName}/fastqc/", mode: 'copy'
    
    container 'tbc'

    input: 
        path readsPath

    output:
        path "${readsPath.baseName}_fastqc.html"
        path "${readsPath.baseName}_fastqc.zip"
    
    script:
    """
    fastqc $readsPath
    """
}